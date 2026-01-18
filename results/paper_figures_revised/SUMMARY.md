# Revised Figure Summary for seqproc Paper

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
