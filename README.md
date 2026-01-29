# seqproc Paper Analysis

Benchmarking and accuracy analysis scripts for the seqproc paper, comparing seqproc against matchbox and splitcode on SPLiT-seq and 10x Genomics data.

## Final Results

The final results used in the paper are located in:
- **`results/paper_figures_revised/`**: Main SPLiT-seq benchmark figures (Runtime, Recovery, Correlation).
- **`results/precision_recall/`**: Precision-recall analysis (Figure E style).
- **`docs/FINAL_REPORT.md`**: Comprehensive summary of the findings.

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Set tool paths (if not in PATH)
export SEQPROC_BIN=/path/to/seqproc
export MATCHBOX_BIN=/path/to/matchbox
export SPLITCODE_BIN=/path/to/splitcode

# 1. Run Unified SPLiT-seq Benchmark (Figures 1, 2, 4)
python scripts/run_splitseq_pipeline.py \
    --r1 /path/to/SRR6750041_1M_R1.fastq \
    --r2 /path/to/SRR6750041_1M_R2.fastq \
    --name SRR6750041_1M \
    --threads 4 \
    --replicates 5

# 2. Run Precision-Recall Analysis (Figure 3)
python scripts/run_precision_recall.py --num-reads 50000 --threads 4

# 3. Run 10x Chromium v2 Benchmark (Splitcode)
python scripts/benchmark_splitcode_10x.py
```

## Scripts

### `run_splitseq_pipeline.py`
The main orchestration script for the SPLiT-seq benchmark. It runs seqproc, matchbox, and splitcode, measures runtime (with replicates), parses output to calculate recovery rates, and compares barcode assignments for correlation analysis.
- **Inputs:** Paired-end FASTQ files.
- **Outputs:** `pipeline_summary.json` and figures in `results/pipeline/<name>/`.

### `run_precision_recall.py`
Generates synthetic data with known ground truth to evaluate precision and recall at different error tolerances.
- **Outputs:** `fig_precision_recall.png` in `results/precision_recall/`.

### `benchmark_splitcode_10x.py`
Specific benchmark for 10x Chromium v2 data using Splitcode with positional extraction, as used for the final paper table updates.

### `generate_revised_figures.py`
Used to generate the polished figures found in `results/paper_figures_revised/` from the raw JSON results.

## Configurations

The analysis relies on specific configuration files for each tool:

### seqproc (`configs/seqproc/`)
- `splitseq_real.geom`: Main SPLiT-seq geometry (Anchor Relative).
- `splitseq_filter.geom`: "Locked-in" geometry with whitelist filtering (referenced in Final Report).
- `10x_longread.geom`: Geometry for 10x Long Read (GridION/PromethION).

### matchbox (`configs/matchbox/`)
- `splitseq.mb`: SPLiT-seq configuration script.
- `10x_longread.mb`: 10x Long Read configuration.

### splitcode (`configs/splitcode/`)
- `splitseq_paper.config`: Main SPLiT-seq configuration.
- `10x_v2_user.config`: Optimized positional extraction config for 10x v2.

## Directory Structure

```
├── configs/                    # Tool configurations (.geom, .mb, .config)
├── data/                       # Input datasets (not included in repo)
├── docs/                       # Final reports and analysis summaries
├── paper_summaries/            # detailed daily summaries of analysis steps
├── results/
│   ├── paper_figures_revised/  # FINAL figures for the paper
│   └── precision_recall/       # Precision-recall analysis results
└── scripts/                    # Analysis and plotting scripts
```

## Citation

If you use this analysis, please cite the seqproc paper.
