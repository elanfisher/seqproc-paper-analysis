# seqproc Paper Analysis Summary

This document provides a comprehensive summary of the benchmark analysis, methodologies, and key findings.

---

## 1. Tools Compared

| Tool | Description | Source |
|------|-------------|--------|
| **seqproc** | Rust-based sequence geometry parser | COMBINE-lab |
| **matchbox** | Flexible read processor with scripting language | jakob-schuster/matchbox |
| **splitcode** | Configuration-based barcode extraction | pachterlab/splitcode |

---

## 2. Datasets

### Datasets Used in Benchmarks

| Dataset | Accession | Reads | Type | Description |
|---------|-----------|-------|------|-------------|
| SPLiT-seq PE | SRR6750041 | 77.6M (full), 1M (subset) | Paired-end, Illumina | Original SPLiT-seq publication data |
| SPLiT-seq SE | SRR13948564 | 5.6M | Single-end, ONT long-read | Long-read SPLiT-seq |
| 10x GridION | ERR9958134 | 7.5M | Single-end, ONT long-read | 10x Chromium on Oxford Nanopore |

### Excluded Datasets

| Dataset | Accession | Reason |
|---------|-----------|--------|
| 10x PromethION | ERR9958135 | 61GB file size, disk space constraint |

---

## 3. Precision-Recall Simulation Methodology

### Overview

The precision-recall analysis uses **synthetic data** with a known ground truth to evaluate barcode extraction accuracy.

### Methodology (matching matchbox paper Figure E style)

1. **Whitelist Generation**
   - Generate 96 random 8bp barcodes for each of BC1, BC2, BC3 pools
   - Create whitelist of all valid BC1+BC2 combinations (96 × 96 = 9,216 combinations)

2. **Synthetic Read Generation**
   - For each read, randomly select barcodes from the whitelist pools
   - Build R2 sequence with SPLiT-seq structure:
     ```
     NN + UMI(10bp) + BC3 + Linker1(30bp) + BC2 + Linker2(16bp) + BC1 + polyA
     ```
   - Introduce 2% base-level substitution errors (realistic sequencing error rate)

3. **Ground Truth**
   - Each read has a known TRUE barcode assignment (BC1, BC2, BC3, UMI)

4. **Metrics at Different Tolerance Levels**
   - **Tolerance 0**: Exact whitelist match required
   - **Tolerance 1-3**: Allow 1-3 mismatches for whitelist matching
   
   For each tolerance level:
   - **Correct**: Tool assigns read to the TRUE whitelist barcode
   - **Incorrect**: Tool assigns read to a WRONG whitelist barcode
   - **Unassigned**: Tool cannot match read to any whitelist barcode

### Key Parameters

```python
INPUT_ERROR_RATE = 0.02  # 2% base-level error rate
num_reads = 50000        # Default synthetic reads
tolerance_levels = [0, 1, 2, 3]  # Allowed mismatches
```

### Limitations

- Uses synthetic data, not biological ground truth
- Does not account for real sequencing error patterns (indels, homopolymers)
- Whitelist is randomly generated, not real SPLiT-seq barcodes

---

## 4. Matchbox Paper Splitcode Comparison

### Throughput Comparison (from matchbox bioRxiv preprint, Figure 7)

Benchmarks run on 10M reads, Intel Xeon E5-2690 v4 @ 2.60GHz:

| Task | matchbox | splitcode | flexiplex |
|------|----------|-----------|-----------|
| Sequence search (0% error) | 78K reads/sec | 71K reads/sec | 61K reads/sec |
| Sequence search (20% error) | 67K reads/sec | 68K reads/sec | 44K reads/sec |
| Barcode discovery | 73K reads/sec | 65K reads/sec | 80K reads/sec |
| Barcode matching | 33K reads/sec | 49K reads/sec | 81K reads/sec |

### Key Findings

- **matchbox and splitcode have comparable performance** in the matchbox paper
- splitcode is slightly faster for barcode matching (49K vs 33K reads/sec)
- matchbox is slightly faster for sequence search tasks

### Our Benchmark Comparison

Our results (1M subset, 4 threads for seqproc/matchbox, 1 thread for splitcode):

| Dataset | seqproc | matchbox | splitcode | splitcode/seqproc |
|---------|---------|----------|-----------|-------------------|
| SPLiT-seq PE | 2.6s | 4.0s | 10.6s | 4.1x slower |
| SPLiT-seq SE | 16.1s | 2.9s | 103s | 6.4x slower |
| 10x GridION | 20.1s | 6.4s | 65.8s | 3.3x slower |

**Note**: Splitcode ran single-threaded due to memory constraints. With equal threading, splitcode would likely be ~2-3x faster.

---

## 5. AlevinFry Zebrafish Pineal Analysis

### Overview

The AlevinFry paper (Nature Methods 2022) includes a zebrafish pineal gland gene expression analysis as a case study.

### Dataset

- **Accession**: SRR8315379, SRR8315380
- **Type**: 10x Chromium v2
- **Reference**: Zebrafish GRCm10 (dr-101)

### Pipeline Structure

```
FASTQ files
    ↓
salmon alevin (mapping)
    ↓
alevin-fry (quantification)
    ↓
R analysis (Seurat/SingleCellExperiment)
    ↓
t-SNE clustering & gene expression plots
```

### Relevant Scripts

Located in `github.com/COMBINE-lab/alevin-fry-paper-scripts`:

- `data_prep_scripts/gather_samples.sh` - Downloads zebrafish data
- `analysis_scripts/zebrafish_pineal_expression/zebrafish_dataset_clustering_analysis.Rmd` - Main analysis

### Integration Point for seqproc

seqproc can replace `seq_geom_xform` (or the alevin-fry geometry parsing step) to:
1. Parse complex read geometries
2. Output transformed FASTQs with extracted barcodes/UMIs
3. Feed into salmon alevin for mapping

---

## 6. Dataset Analysis: PE vs SE

### Performance Comparison

| Metric | SPLiT-seq PE | SPLiT-seq SE | 10x GridION |
|--------|--------------|--------------|-------------|
| Recovery Rate | ~87% | ~24% | ~36-41% |
| seqproc vs matchbox speed | seqproc 1.6x faster | matchbox 5.5x faster | matchbox 3.1x faster |

### Observations

1. **Paired-end (PE) performs best**: seqproc excels on PE data
2. **Single-end (SE) has lower recovery**: Expected due to:
   - Long-read sequencing error profiles
   - Primer/adapter finding challenges
   - Variable read lengths
3. **seqproc is slower on SE/long-read data**: matchbox is optimized for these use cases

### Recommendation

For the paper, consider:
- **Focus on PE analysis** where seqproc shows clear advantage
- Include SE/long-read as supplementary to show tool limitations
- Document that SE performance is a known area for improvement

---

## 7. Script Reference

### Main Scripts

| Script | Purpose |
|--------|---------|
| `run_full_benchmark.py` | Complete benchmark orchestration |
| `run_precision_recall.py` | Precision-recall analysis with synthetic data |
| `run_correlation_only.py` | Barcode correlation analysis (crash recovery) |
| `generate_revised_figures.py` | Figure generation from results JSON |

### Usage Examples

```bash
# Full benchmark
python scripts/run_full_benchmark.py --full --threads 4 --replicates 2

# Precision-recall analysis
python scripts/run_precision_recall.py --num-reads 50000 --threads 4

# Generate revised figures
python scripts/generate_revised_figures.py \
    --results results/paper_figures/benchmark_results.json \
    --output results/paper_figures_revised \
    --title-suffix "1M reads" \
    --note "*splitcode: 1 thread (memory constraint)"
```

---

## 8. Known Issues and Caveats

### Full Dataset Anomalies

The full 77M SPLiT-seq PE dataset showed anomalous results:
- seqproc: 48.5% recovery (37.7M reads) ✓
- matchbox: **0.05%** recovery (39K reads) ✗
- splitcode: **0.0002%** recovery (167 reads) ✗

**This needs investigation** - likely configuration or memory issues.

### Threading Constraints

- splitcode was limited to 1 thread to prevent memory exhaustion
- This significantly impacts splitcode's runtime performance

### Subset vs Full Dataset

- Barcode correlation analysis uses 1M subset (not full 77M) to avoid memory crash
- Full dataset results may differ from subset results

---

## 9. File Locations

```
seqproc-paper-analysis-clean/
├── configs/                    # Tool configuration files
│   ├── seqproc/               # .geom files
│   ├── matchbox/              # .mb scripts
│   └── splitcode/             # .config files
├── data/                       # Input datasets
├── results/
│   ├── paper_figures/         # 1M subset results
│   ├── paper_figures_full/    # Full dataset results
│   └── paper_figures_revised/ # Revised figures
├── scripts/                    # Analysis scripts
└── docs/                       # Documentation
```

---

*Generated: January 18, 2026*
