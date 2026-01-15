# seqproc Paper Analysis

Benchmarking and accuracy analysis scripts for the seqproc paper, comparing seqproc against matchbox and splitcode on SPLiT-seq data.

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run unified pipeline (benchmark + accuracy)
python scripts/run_splitseq_pipeline.py \
    --r1 /path/to/R1.fastq \
    --r2 /path/to/R2.fastq \
    --name my_dataset \
    --threads 4 \
    --replicates 5
```

## Results

The pipeline generates:
- `fig_runtime_boxplot.png` - Runtime comparison box plot
- `fig_speedup.png` - Speedup vs matchbox baseline
- `fig_match_rate.png` - Read match rate comparison
- `fig_correlation.png` - Barcode count correlation (R²)
- `pipeline_summary.json` - All metrics

## Tool Dependencies

Install the following tools and update paths in `scripts/run_splitseq_pipeline.py`:

- [seqproc](https://github.com/COMBINE-lab/seqproc) - Fast sequence preprocessor
- [matchbox](https://github.com/noamteyssier/matchbox) - Pattern matching tool
- [splitcode](https://github.com/pachterlab/splitcode) - Barcode demultiplexer

## Data

SPLiT-seq benchmark data from SRA:
- `SRR6750041` - 77.6M paired-end reads (SPLiT-seq Round 2)

## Figures

Matches the style of figures from the matchbox paper:
- **Runtime box plots** with multiple replicates
- **Barcode count correlation** (Fig H style, R² metric)
- **Read match rate** comparison across tools

## Directory Structure

```
├── scripts/
│   └── run_splitseq_pipeline.py  # Main unified pipeline
├── configs/
│   ├── seqproc/                  # seqproc geometry files
│   ├── matchbox/                 # matchbox scripts (reference)
│   └── splitcode/                # splitcode config files
└── results/                      # Output directory
```

## Citation

If you use this analysis, please cite the seqproc paper.
