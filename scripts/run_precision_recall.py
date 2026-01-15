#!/usr/bin/env python3
"""
Precision-Recall Analysis (Matchbox Paper Fig E style)

Generates synthetic SPLiT-seq data with known barcodes, then measures
precision and recall at different error rates for each tool.

Usage:
    python scripts/run_precision_recall.py --num-reads 50000 --threads 4
"""

import os
import sys
import subprocess
import tempfile
import json
import argparse
import random
from pathlib import Path
from collections import defaultdict
from typing import Dict, Tuple, List
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).parent.parent
os.chdir(PROJECT_ROOT)

# Tool binaries
SEQPROC_BIN = os.environ.get("SEQPROC_BIN", "seqproc")
MATCHBOX_BIN = os.environ.get("MATCHBOX_BIN", "matchbox")
SPLITCODE_BIN = os.environ.get("SPLITCODE_BIN", "splitcode")

# Config files
SEQPROC_CONFIG = "configs/seqproc/splitseq_real.geom"
MATCHBOX_CONFIG = "configs/matchbox/splitseq.mb"
SPLITCODE_CONFIG = "configs/splitcode/splitseq_paper.config"

# SPLiT-seq structure
LINKER1 = "GTGGCCGCTGTTTCGCATCGGCGTACGACT"  # 30bp
LINKER2 = "ATCCACGTGCTTGAGA"  # 16bp

COLORS = {
    'seqproc': '#4DBBD5',
    'matchbox': '#E64B35',
    'splitcode': '#3C5488',
}

MARKERS = {0: 'o', 1: 's', 2: '^', 3: 'D'}


def random_barcode(length: int = 8) -> str:
    """Generate random barcode sequence."""
    return ''.join(random.choice('ACGT') for _ in range(length))


def introduce_errors(seq: str, error_rate: float) -> str:
    """Introduce random substitution errors."""
    if error_rate <= 0:
        return seq
    result = list(seq)
    for i in range(len(result)):
        if random.random() < error_rate:
            bases = [b for b in 'ACGT' if b != result[i]]
            result[i] = random.choice(bases)
    return ''.join(result)


def generate_synthetic_data(num_reads: int, error_rate: float, output_dir: str) -> Dict[str, Tuple]:
    """Generate synthetic SPLiT-seq data with known barcodes.
    
    Returns dict mapping read_id -> (bc1, bc2, bc3, umi)
    """
    random.seed(42)  # Reproducibility
    
    # Generate barcode pools
    bc1_pool = [random_barcode(8) for _ in range(96)]
    bc2_pool = [random_barcode(8) for _ in range(96)]
    bc3_pool = [random_barcode(8) for _ in range(96)]
    
    ground_truth = {}
    
    r1_path = f"{output_dir}/synthetic_R1.fq"
    r2_path = f"{output_dir}/synthetic_R2.fq"
    
    with open(r1_path, 'w') as r1_out, open(r2_path, 'w') as r2_out:
        for i in range(num_reads):
            read_id = f"SYN.{i+1}"
            
            # Select random barcodes
            bc1 = random.choice(bc1_pool)
            bc2 = random.choice(bc2_pool)
            bc3 = random.choice(bc3_pool)
            umi = random_barcode(10)
            
            ground_truth[read_id] = (bc1, bc2, bc3, umi)
            
            # Build R2 sequence: NN + UMI + BC3 + Linker1 + BC2 + Linker2 + BC1 + polyA
            r2_seq = (
                "NN" +  # Skip bases
                introduce_errors(umi, error_rate) +
                introduce_errors(bc3, error_rate) +
                introduce_errors(LINKER1, error_rate) +
                introduce_errors(bc2, error_rate) +
                introduce_errors(LINKER2, error_rate) +
                introduce_errors(bc1, error_rate) +
                "AAAAAAAAAA"  # polyA tail
            )
            
            # R1 is cDNA (random sequence)
            r1_seq = ''.join(random.choice('ACGT') for _ in range(100))
            
            # Quality scores (all high quality)
            r1_qual = 'I' * len(r1_seq)
            r2_qual = 'I' * len(r2_seq)
            
            # Write FASTQ
            r1_out.write(f"@{read_id}\n{r1_seq}\n+\n{r1_qual}\n")
            r2_out.write(f"@{read_id}\n{r2_seq}\n+\n{r2_qual}\n")
    
    return ground_truth, r1_path, r2_path


def run_seqproc(r1_path: str, r2_path: str, threads: int, tmpdir: str) -> Dict[str, Tuple]:
    """Run seqproc and extract barcodes."""
    out1 = f"{tmpdir}/seqproc_R1.fq"
    out2 = f"{tmpdir}/seqproc_R2.fq"
    cmd = f"{SEQPROC_BIN} --geom {SEQPROC_CONFIG} --file1 {r1_path} --file2 {r2_path} --out1 {out1} --out2 {out2} --threads {threads}"
    subprocess.run(cmd, shell=True, capture_output=True)
    
    results = {}
    if os.path.exists(out2):
        with open(out2, 'r') as f:
            while True:
                header = f.readline()
                if not header:
                    break
                seq = f.readline().strip()
                f.readline()
                f.readline()
                
                read_id = header.strip().split()[0].replace('@', '')
                if len(seq) >= 32:
                    umi = seq[0:10]
                    bc3 = seq[10:18]
                    bc1 = seq[-8:]
                    bc2 = seq[-16:-8]
                    results[read_id] = (bc1, bc2, bc3, umi)
    return results


def run_matchbox(r2_path: str, threads: int, tmpdir: str) -> Dict[str, Tuple]:
    """Run matchbox and extract barcodes."""
    out_tsv = f"{tmpdir}/matchbox_out.tsv"
    with open(MATCHBOX_CONFIG, 'r') as f:
        mb_script = f.read()
    cmd = f'{MATCHBOX_BIN} -e 0.2 -t {threads} -r "{mb_script}" {r2_path} > {out_tsv}'
    subprocess.run(cmd, shell=True, capture_output=True)
    
    results = {}
    if os.path.exists(out_tsv):
        with open(out_tsv, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    read_id = parts[0]
                    bc1, bc2, bc3, umi = parts[1], parts[2], parts[3], parts[4]
                    results[read_id] = (bc1, bc2, bc3, umi)
    return results


def run_splitcode(r1_path: str, r2_path: str, threads: int, tmpdir: str) -> Dict[str, Tuple]:
    """Run splitcode and extract barcodes."""
    out1 = f"{tmpdir}/splitcode_R1.fq"
    out2 = f"{tmpdir}/splitcode_R2.fq"
    mapping = f"{tmpdir}/splitcode_mapping.txt"
    cmd = f"{SPLITCODE_BIN} -c {SPLITCODE_CONFIG} --assign -N 2 -t {threads} -m {mapping} -o {out1},{out2} {r1_path} {r2_path}"
    subprocess.run(cmd, shell=True, capture_output=True)
    
    # Splitcode doesn't output sequences, just assignments - count matched reads
    results = {}
    if os.path.exists(out2):
        with open(out2, 'r') as f:
            while True:
                header = f.readline()
                if not header:
                    break
                seq = f.readline().strip()
                f.readline()
                f.readline()
                read_id = header.strip().split()[0].replace('@', '')
                # Mark as matched (we can't extract actual barcodes from splitcode)
                results[read_id] = ('', '', '', '')
    return results


def calculate_precision_recall(ground_truth: Dict, predicted: Dict, allow_mismatch: int = 1) -> Tuple[float, float]:
    """Calculate precision and recall.
    
    - Correct: predicted barcode matches ground truth (within allow_mismatch)
    - Incorrect: predicted barcode doesn't match ground truth
    - Missed: ground truth not in predicted
    """
    correct = 0
    incorrect = 0
    
    for read_id, pred_bc in predicted.items():
        if read_id in ground_truth:
            true_bc = ground_truth[read_id]
            # Compare BC1, BC2, BC3 (not UMI since it varies)
            matches = True
            for i in range(3):  # bc1, bc2, bc3
                if pred_bc[i] and true_bc[i]:
                    dist = sum(c1 != c2 for c1, c2 in zip(pred_bc[i], true_bc[i]))
                    if dist > allow_mismatch:
                        matches = False
                        break
            if matches:
                correct += 1
            else:
                incorrect += 1
        else:
            incorrect += 1
    
    total_predicted = len(predicted)
    total_truth = len(ground_truth)
    
    # Fraction of reads with correct barcodes (of those predicted)
    frac_correct = correct / total_truth if total_truth > 0 else 0
    # Fraction of reads with incorrect barcodes (of those predicted)
    frac_incorrect = incorrect / total_truth if total_truth > 0 else 0
    
    return frac_correct, frac_incorrect


def main():
    parser = argparse.ArgumentParser(description="Precision-Recall Analysis")
    parser.add_argument('--num-reads', type=int, default=50000, help='Number of synthetic reads')
    parser.add_argument('--threads', type=int, default=4, help='Threads per tool')
    parser.add_argument('--output-dir', default=None)
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir) if args.output_dir else PROJECT_ROOT / "results/precision_recall"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("Precision-Recall Analysis (Fig E style)")
    print("=" * 70)
    print(f"Synthetic reads: {args.num_reads:,}")
    print(f"Threads: {args.threads}")
    print(f"Output: {output_dir}")
    print("=" * 70)
    
    # Error rates to test (0, 1, 2, 3 errors per 100bp ~ 0%, 1%, 2%, 3%)
    error_rates = [0.0, 0.01, 0.02, 0.03]
    tools = ['seqproc', 'matchbox', 'splitcode']
    
    results = {tool: {'frac_correct': [], 'frac_incorrect': []} for tool in tools}
    
    for error_idx, error_rate in enumerate(error_rates):
        print(f"\n[{error_idx+1}/{len(error_rates)}] Error rate: {error_rate:.0%}")
        
        with tempfile.TemporaryDirectory() as tmpdir:
            # Generate synthetic data
            print("  Generating synthetic data...", end=" ", flush=True)
            ground_truth, r1_path, r2_path = generate_synthetic_data(
                args.num_reads, error_rate, tmpdir
            )
            print(f"{len(ground_truth):,} reads")
            
            # Run seqproc
            print("  Running seqproc...", end=" ", flush=True)
            sp_pred = run_seqproc(r1_path, r2_path, args.threads, tmpdir)
            sp_correct, sp_incorrect = calculate_precision_recall(ground_truth, sp_pred)
            results['seqproc']['frac_correct'].append(sp_correct)
            results['seqproc']['frac_incorrect'].append(sp_incorrect)
            print(f"correct={sp_correct:.3f}, incorrect={sp_incorrect:.3f}")
            
            # Run matchbox
            print("  Running matchbox...", end=" ", flush=True)
            mb_pred = run_matchbox(r2_path, args.threads, tmpdir)
            mb_correct, mb_incorrect = calculate_precision_recall(ground_truth, mb_pred)
            results['matchbox']['frac_correct'].append(mb_correct)
            results['matchbox']['frac_incorrect'].append(mb_incorrect)
            print(f"correct={mb_correct:.3f}, incorrect={mb_incorrect:.3f}")
            
            # Run splitcode (can only measure recall, not barcode accuracy)
            print("  Running splitcode...", end=" ", flush=True)
            sc_pred = run_splitcode(r1_path, r2_path, args.threads, tmpdir)
            # For splitcode, we measure read match rate only
            sc_match_rate = len(sc_pred) / len(ground_truth)
            results['splitcode']['frac_correct'].append(sc_match_rate)
            results['splitcode']['frac_incorrect'].append(1 - sc_match_rate)
            print(f"match_rate={sc_match_rate:.3f}")
    
    # Generate Fig E style plot
    print("\nGenerating precision-recall plot...")
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    for tool in ['seqproc', 'matchbox']:  # Only tools with barcode-level accuracy
        frac_correct = results[tool]['frac_correct']
        frac_incorrect = results[tool]['frac_incorrect']
        
        # Plot line
        ax.plot(frac_incorrect, frac_correct, 
                color=COLORS[tool], linewidth=2, label=tool, alpha=0.8)
        
        # Plot points with different markers for error rates
        for i, (fi, fc) in enumerate(zip(frac_incorrect, frac_correct)):
            ax.scatter(fi, fc, color=COLORS[tool], marker=MARKERS[i], 
                      s=100, edgecolor='black', linewidth=1, zorder=5)
    
    # Add legend for error rates
    for i, er in enumerate(error_rates):
        ax.scatter([], [], color='gray', marker=MARKERS[i], s=80, 
                  label=f'Error rate {int(er*100)}', edgecolor='black')
    
    ax.set_xlabel('Fraction of reads with incorrect barcodes', fontsize=12)
    ax.set_ylabel('Fraction of reads with correct barcodes', fontsize=12)
    ax.set_title('Precision-Recall Analysis\n(Synthetic SPLiT-seq Data)', fontsize=14, fontweight='bold')
    ax.legend(loc='lower right', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-0.005, max(0.1, max(max(results[t]['frac_incorrect']) for t in tools) * 1.1))
    ax.set_ylim(0, 1.05)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'fig_precision_recall.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'fig_precision_recall.pdf', bbox_inches='tight')
    plt.close()
    
    # Save results
    with open(output_dir / 'precision_recall_results.json', 'w') as f:
        json.dump({
            'error_rates': error_rates,
            'results': results,
            'num_reads': args.num_reads,
        }, f, indent=2)
    
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"\n{'Error Rate':<12} {'seqproc':>20} {'matchbox':>20}")
    print("-" * 54)
    for i, er in enumerate(error_rates):
        sp_c = results['seqproc']['frac_correct'][i]
        mb_c = results['matchbox']['frac_correct'][i]
        print(f"{er*100:>3.0f}%         {sp_c*100:>18.1f}% {mb_c*100:>18.1f}%")
    print("-" * 54)
    print(f"\nFigure saved to: {output_dir / 'fig_precision_recall.png'}")
    print("=" * 70)


if __name__ == "__main__":
    main()
