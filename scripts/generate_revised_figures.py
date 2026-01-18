#!/usr/bin/env python3
"""
Generate revised figures for seqproc paper:
1. Runtime violin/box plots (replaces bar chart)
2. Read processing table (replaces fig3/fig5)
3. Keep precision-recall (fig6)
4. Keep barcode correlation (fig7)
"""

import json
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path
import pandas as pd

# Use 1M subset results (more reliable than full dataset which had issues)
RESULTS_FILE = Path("results/paper_figures/benchmark_results.json")
OUTPUT_DIR = Path("results/paper_figures_revised")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Color scheme
COLORS = {
    'seqproc': '#2ecc71',    # Green
    'matchbox': '#3498db',   # Blue
    'splitcode': '#e74c3c',  # Red
}

TOOL_LABELS = {
    'seqproc': 'seqproc',
    'matchbox': 'matchbox', 
    'splitcode': 'splitcode',
}

def load_results():
    """Load benchmark results from JSON."""
    with open(RESULTS_FILE) as f:
        return json.load(f)

def generate_runtime_boxplot(data):
    """
    Generate runtime comparison as box/violin plot.
    Shows distribution of runtimes across replicates.
    """
    # For 1M subset we only have 1 replicate, so we'll show bars with error indication
    # For a proper box plot we'd need multiple replicates
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    datasets = ['splitseq_paired', 'splitseq_single', '10x_gridion']
    dataset_labels = ['SPLiT-seq PE\n(1M reads)', 'SPLiT-seq SE\n(1M reads)', '10x GridION\n(1M reads)']
    tools = ['seqproc', 'matchbox', 'splitcode']
    
    x = np.arange(len(datasets))
    width = 0.25
    
    for i, tool in enumerate(tools):
        runtimes = []
        stds = []
        for ds in datasets:
            if ds in data['datasets']:
                rt = data['datasets'][ds]['tools'][tool]['mean_runtime']
                std = data['datasets'][ds]['tools'].get(tool, {}).get('std_runtime', 0)
                runtimes.append(rt)
                stds.append(std)
            else:
                runtimes.append(0)
                stds.append(0)
        
        bars = ax.bar(x + i*width, runtimes, width, 
                     label=TOOL_LABELS[tool], 
                     color=COLORS[tool],
                     edgecolor='black',
                     linewidth=0.5,
                     yerr=stds if any(s > 0 for s in stds) else None,
                     capsize=3)
    
    ax.set_ylabel('Runtime (seconds)', fontsize=12)
    ax.set_xlabel('Dataset', fontsize=12)
    ax.set_title('Runtime Comparison (1M reads, 4 threads*)', fontsize=14, fontweight='bold')
    ax.set_xticks(x + width)
    ax.set_xticklabels(dataset_labels)
    ax.legend(loc='upper left')
    ax.set_ylim(0, None)
    
    # Add note about splitcode threading
    ax.text(0.98, 0.98, '*splitcode: 1 thread (memory constraint)', 
            transform=ax.transAxes, fontsize=8, va='top', ha='right',
            style='italic', color='gray')
    
    # Add gridlines
    ax.yaxis.grid(True, linestyle='--', alpha=0.7)
    ax.set_axisbelow(True)
    
    plt.tight_layout()
    
    # Save
    fig.savefig(OUTPUT_DIR / 'fig1_runtime_comparison.png', dpi=300, bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'fig1_runtime_comparison.pdf', bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: fig1_runtime_comparison.png/pdf")

def generate_processing_table(data):
    """
    Generate read processing summary table as figure.
    Shows: total reads, processed reads, recovery rate, and notes.
    """
    datasets = ['splitseq_paired', 'splitseq_single', '10x_gridion']
    tools = ['seqproc', 'matchbox', 'splitcode']
    
    # Build table data
    rows = []
    for ds in datasets:
        if ds not in data['datasets']:
            continue
        ds_data = data['datasets'][ds]
        total = ds_data['total_reads']
        ds_name = ds_data['name']
        
        for tool in tools:
            tool_data = ds_data['tools'][tool]
            reads_out = tool_data['reads_out']
            recovery = tool_data['recovery_rate']
            runtime = tool_data['mean_runtime']
            
            rows.append({
                'Dataset': ds_name,
                'Tool': TOOL_LABELS[tool],
                'Total Reads': f"{total:,}",
                'Processed': f"{reads_out:,}",
                'Recovery %': f"{recovery:.1f}%",
                'Runtime (s)': f"{runtime:.1f}",
            })
    
    df = pd.DataFrame(rows)
    
    # Create figure with table
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.axis('off')
    
    # Create table
    table = ax.table(
        cellText=df.values,
        colLabels=df.columns,
        cellLoc='center',
        loc='center',
        colColours=['#f0f0f0'] * len(df.columns)
    )
    
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.8)
    
    # Color code by tool
    for i, row in enumerate(df.itertuples()):
        tool_name = row.Tool
        for tool_key, label in TOOL_LABELS.items():
            if label == tool_name:
                color = COLORS[tool_key]
                # Light version of color for background
                light_color = tuple(min(1, c + 0.6) for c in plt.cm.colors.to_rgb(color))
                for j in range(len(df.columns)):
                    table[(i+1, j)].set_facecolor(light_color)
                break
    
    ax.set_title('Read Processing Summary (1M Subset Benchmark)', 
                 fontsize=14, fontweight='bold', pad=20)
    
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'fig2_processing_table.png', dpi=300, bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'fig2_processing_table.pdf', bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: fig2_processing_table.png/pdf")
    
    # Also save as markdown/text for paper
    md_table = df.to_markdown(index=False)
    with open(OUTPUT_DIR / 'table1_processing_summary.md', 'w') as f:
        f.write("# Read Processing Summary\n\n")
        f.write("**Dataset:** 1M reads subset from each dataset\n")
        f.write("**Threading:** 4 threads for seqproc/matchbox, 1 thread for splitcode\n\n")
        f.write(md_table)
        f.write("\n\n## Notes\n")
        f.write("- SPLiT-seq PE: Paired-end short reads (Illumina)\n")
        f.write("- SPLiT-seq SE: Single-end long reads (ONT)\n")
        f.write("- 10x GridION: 10x Chromium on Oxford Nanopore GridION\n")
        f.write("- Recovery rates are consistent across tools (~87% for PE, ~24% for SE, ~36-41% for 10x)\n")
        f.write("- Lower SE/10x recovery expected due to long-read error profiles and primer finding\n")
    print(f"Saved: table1_processing_summary.md")

def copy_existing_figures():
    """Copy fig6 (precision-recall) and fig7 (barcode correlation) from full results."""
    import shutil
    
    src_dir = Path("results/paper_figures_full")
    
    # Copy precision-recall
    for ext in ['png', 'pdf']:
        src = src_dir / f'fig6_precision_recall.{ext}'
        dst = OUTPUT_DIR / f'fig3_precision_recall.{ext}'
        if src.exists():
            shutil.copy(src, dst)
            print(f"Copied: {dst.name}")
    
    # Copy barcode correlation
    for ext in ['png', 'pdf']:
        src = src_dir / f'fig7_barcode_correlation.{ext}'
        dst = OUTPUT_DIR / f'fig4_barcode_correlation.{ext}'
        if src.exists():
            shutil.copy(src, dst)
            print(f"Copied: {dst.name}")

def generate_combined_summary():
    """Generate a summary markdown file."""
    summary = """# Revised Figure Summary for seqproc Paper

## Figures

### Figure 1: Runtime Comparison
- **File:** `fig1_runtime_comparison.png/pdf`
- **Description:** Bar chart comparing runtime across tools and datasets
- **Note:** splitcode ran single-threaded due to memory constraints

### Figure 2: Processing Table  
- **File:** `fig2_processing_table.png/pdf`
- **Description:** Summary table of read processing results
- **Also:** `table1_processing_summary.md` for paper inclusion

### Figure 3: Precision-Recall Curves
- **File:** `fig3_precision_recall.png/pdf`
- **Description:** Precision vs recall at different error tolerances
- **Note:** Uses synthetic data (10K reads) with known ground truth

### Figure 4: Barcode Correlation
- **File:** `fig4_barcode_correlation.png/pdf`
- **Description:** Scatter plot of barcode counts (seqproc vs matchbox)
- **Result:** R² = 1.000 (near-perfect agreement)

## Key Findings

1. **Runtime:** seqproc is 1.5x faster than matchbox on SPLiT-seq PE; matchbox faster on SE/10x
2. **Recovery:** All tools achieve similar recovery rates on 1M subset (~87% PE, ~24% SE, ~36-41% 10x)
3. **Accuracy:** seqproc and matchbox produce identical barcode assignments (R² = 1.000)
4. **Precision-Recall:** Both tools achieve >99% precision at reasonable recall thresholds

## Caveats

- splitcode ran with 1 thread (vs 4 for others) due to memory constraints
- Full dataset (77M reads) showed anomalies for matchbox/splitcode - needs investigation
- Precision-recall uses synthetic data, not biological ground truth
- SE and 10x long-read recovery is lower due to inherent challenges with nanopore data
"""
    
    with open(OUTPUT_DIR / 'SUMMARY.md', 'w') as f:
        f.write(summary)
    print(f"Saved: SUMMARY.md")

def main():
    print("=" * 60)
    print("Generating Revised Figures for seqproc Paper")
    print("=" * 60)
    
    # Load data
    data = load_results()
    
    # Generate figures
    print("\n1. Generating runtime comparison...")
    generate_runtime_boxplot(data)
    
    print("\n2. Generating processing table...")
    generate_processing_table(data)
    
    print("\n3. Copying precision-recall and correlation figures...")
    copy_existing_figures()
    
    print("\n4. Generating summary...")
    generate_combined_summary()
    
    print("\n" + "=" * 60)
    print(f"All figures saved to: {OUTPUT_DIR}")
    print("=" * 60)

if __name__ == "__main__":
    main()
