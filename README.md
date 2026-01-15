# seqproc Paper Analysis

Benchmarking and accuracy analysis scripts for the seqproc paper, comparing seqproc against matchbox and splitcode on SPLiT-seq data.

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Set tool paths
export SEQPROC_BIN=/path/to/seqproc
export MATCHBOX_BIN=/path/to/matchbox
export SPLITCODE_BIN=/path/to/splitcode

# Run unified pipeline (benchmark + accuracy)
python scripts/run_splitseq_pipeline.py \
    --r1 /path/to/R1.fastq \
    --r2 /path/to/R2.fastq \
    --name my_dataset \
    --threads 4 \
    --replicates 5

# Run precision-recall analysis (Fig E style)
python scripts/run_precision_recall.py --num-reads 50000 --threads 4
```

## Scripts

### `run_splitseq_pipeline.py` - Unified Benchmark + Accuracy Pipeline
Generates:
- `fig_runtime_boxplot.png` - Runtime comparison box plot with replicates
- `fig_speedup.png` - Speedup vs matchbox baseline
- `fig_match_rate.png` - Read match rate comparison (Fig C style)
- `fig_correlation_seqproc_matchbox.png` - Barcode count correlation (Fig H style)
- `fig_correlation_seqproc_splitcode.png` - Barcode count correlation
- `pipeline_summary.json` - All metrics

### `run_precision_recall.py` - Precision-Recall Analysis
Generates synthetic data with known barcodes and measures accuracy at different error rates (Fig E style).

## Tool Dependencies

Set paths via environment variables:

- `SEQPROC_BIN` - [seqproc](https://github.com/COMBINE-lab/seqproc)
- `MATCHBOX_BIN` - [matchbox](https://github.com/noamteyssier/matchbox)
- `SPLITCODE_BIN` - [splitcode](https://github.com/pachterlab/splitcode)

## Data

SPLiT-seq benchmark data from SRA:
- `SRR6750041` - 77.6M paired-end reads (SPLiT-seq Round 2)

## Figures (Matchbox Paper Style)

| Figure | Script | Description |
|--------|--------|-------------|
| Fig C | `run_splitseq_pipeline.py` | Read match rate bar chart |
| Fig E | `run_precision_recall.py` | Precision-recall at different error rates |
| Fig H | `run_splitseq_pipeline.py` | Barcode count correlation (R²) |

## Directory Structure

```
├── scripts/
│   ├── run_splitseq_pipeline.py  # Main unified pipeline
│   └── run_precision_recall.py   # Precision-recall analysis
├── configs/
│   ├── seqproc/splitseq_real.geom
│   ├── matchbox/splitseq.mb
│   └── splitcode/splitseq_paper.config
└── results/                      # Output directory
```

## Citation

If you use this analysis, please cite the seqproc paper.
