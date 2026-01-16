# Full Dataset Benchmark Results

This directory contains benchmark results comparing **seqproc**, **matchbox**, and **splitcode** on full datasets.

## Datasets Used

| Dataset | Accession | Reads | Description |
|---------|-----------|-------|-------------|
| SPLiT-seq Paired-End | SRR6750041 | 77.6M | Short-read Illumina paired-end |
| SPLiT-seq Single-End | SRR13948564 | 5.6M | Long-read single-end (ONT) |
| 10x GridION | ERR9958134 | 7.5M | 10x Chromium on Oxford Nanopore GridION |

**Note:** 10x PromethION (ERR9958135) was excluded due to disk space constraints (61GB SRA file).

## Generated Figures

| Figure | Description |
|--------|-------------|
| `fig1_runtime_comparison.png` | Runtime comparison bar chart for all tools across datasets |
| `fig2_speedup_chart.png` | seqproc speedup relative to matchbox |
| `fig3_read_recovery.png` | Percentage of reads successfully processed by each tool |
| `fig4_summary_dashboard.png` | 4-panel summary combining key metrics |
| `fig5_match_rate.png` | Match rate per dataset for each tool |
| `fig6_precision_recall.png` | Precision-recall curves at different error tolerances |
| `fig7_barcode_correlation.png` | Barcode count correlation between seqproc and matchbox |

## Benchmark Configuration

- **Threads:** 4 for seqproc/matchbox, 1 for splitcode (memory safety)
- **Replicates:** 2 per dataset
- **Hardware:** AWS EC2 instance with 32GB RAM

---

## Caveats and Assumptions

### 1. SPLiT-seq Paired-End Results (Critical Issue)

The SPLiT-seq paired-end full dataset (77.6M reads) shows **anomalous results** for matchbox and splitcode:

```
seqproc:   ~38M reads recovered (49%)
matchbox:  ~77K reads recovered (0.1%)  ← Unexpectedly low
splitcode: ~315 reads recovered (0%)    ← Unexpectedly low
```

**Possible causes:**
- Memory pressure causing early termination
- Configuration mismatch for full dataset processing
- Tool-specific issues with very large files

**Recommendation:** These results should be verified on a machine with more RAM (64GB+) or by processing in chunks.

### 2. Splitcode Thread Limitation

Splitcode was run with **1 thread only** to prevent memory exhaustion. This significantly impacts its runtime performance:

- Splitcode is ~10x slower than with full threading
- Runtime comparisons are not apples-to-apples for splitcode
- This was necessary to complete the benchmark without crashing

### 3. Precision-Recall Analysis

- Uses **synthetic data** (10,000 reads) generated from the SPLiT-seq barcode whitelist
- Errors are introduced randomly to simulate sequencing errors
- Measures how well tools recover the "true" barcode at different error tolerance levels
- **Does not use real biological data** for this analysis

### 4. Barcode Correlation Analysis

- Uses **1M read subset** (not full dataset) to avoid memory crash
- The correlation (R² = 1.000) demonstrates that seqproc and matchbox produce identical barcode assignments
- Full dataset correlation was not completed due to memory constraints

### 5. Read Recovery Rate Interpretation

Read recovery rates vary significantly by dataset and tool:

- **SPLiT-seq PE:** ~87% recovery (both tools agree well on 1M subset)
- **SPLiT-seq SE:** ~24% recovery (primer-based extraction is stricter)
- **10x GridION:** 36-41% recovery (long-read specific challenges)

Low recovery rates are expected for long-read data due to:
- Higher error rates in nanopore sequencing
- Variable read lengths
- Adapter/primer finding challenges

### 6. Splitcode Configuration Notes

**SPLiT-seq Single-End:** Custom config (`splitseq_singleend.config`) extracts barcodes after primer anchor.

**10x Long-Read:** Splitcode uses a simplified extraction pattern that doesn't perform whitelist matching (unlike seqproc/matchbox which use the 10x barcode whitelist). This means:
- Splitcode extracts ANY 16bp sequence at the expected position
- seqproc/matchbox only output reads matching a valid 10x barcode
- Read counts are comparable but biological validity differs

### 7. Runtime Measurement

- Runtimes include I/O overhead (reading/writing FASTQ files)
- Temporary files are written to `/tmp` (local SSD)
- No warm-up runs were performed

---

## Shortcuts and Corners Cut

1. **Splitcode single-threaded:** Due to memory constraints, splitcode ran with 1 thread vs 4 for other tools.

2. **PromethION dataset excluded:** The 61GB ERR9958135 dataset couldn't fit on disk.

3. **Correlation on subset:** Barcode correlation used 1M reads instead of full 77M to avoid crashes.

4. **2 replicates only:** Used 2 replicates instead of typical 3-5 for statistical power.

5. **No biological validation:** Results show computational agreement, not biological accuracy.

6. **Synthetic PR data:** Precision-recall uses synthetic reads, not validated ground truth.

---

## Reproducing These Results

```bash
# Full benchmark (requires ~68GB disk space and 32GB+ RAM)
python scripts/run_full_benchmark.py --full --threads 4 --splitcode-threads 1 --replicates 2

# Just the correlation analysis (if crashed)
python scripts/run_correlation_only.py --threads 4
```

## File Manifest

```
results/paper_figures_full/
├── README.md                      # This file
├── RESULTS.md                     # Auto-generated summary
├── benchmark_results.json         # Raw benchmark data
├── fig1_runtime_comparison.png    # Runtime comparison
├── fig2_speedup_chart.png         # Speedup analysis
├── fig3_read_recovery.png         # Read recovery rates
├── fig4_summary_dashboard.png     # Summary dashboard
├── fig5_match_rate.png            # Match rates
├── fig6_precision_recall.png      # Precision-recall curves
└── fig7_barcode_correlation.png   # Barcode correlation
```

---

*Generated: January 16, 2026*
