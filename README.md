# seqproc Paper Analysis

Benchmarking and accuracy analysis scripts for the seqproc paper, comparing seqproc against matchbox and splitcode on SPLiT-seq and 10x Genomics data.

## Final Results

The final results used in the paper are located in:
- **`results/paper_figures_revised/`**: Main SPLiT-seq benchmark figures (Runtime, Recovery, Correlation).
- **`results/paper_figures/`**: Comprehensive benchmark tables (Recovery, Performance) across all datasets (SPLiT-seq PE, SE, 10x Long Read).
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
```

## Reproduction Steps

The following scripts were used to generate the data for the final paper draft.

### 1. Main Paper Benchmarks (Tables & Figures)
This script runs the core benchmarks across multiple datasets (SPLiT-seq Paired-End, SPLiT-seq Single-End, 10x GridION, 10x PromethION) to generate the performance distribution and recovery tables.

**Script:** `scripts/run_paper_benchmarks.py`
**Output:** `results/paper_figures/` (Figures 1, 2, 3 and `benchmark_results.json`)

```bash
# Run benchmarks with 3 replicates (requires ~1M reads/dataset in data/)
python scripts/run_paper_benchmarks.py --threads 4 --replicates 3
```

**Datasets involved:**
- SPLiT-seq PE (`SRR6750041`)
- SPLiT-seq SE Long Read (`SRR13948564`)
- 10x GridION (`ERR9958134`)
- 10x PromethION (`ERR9958135`)

### 2. Precision-Recall Analysis (Figure 3/E)
Generates synthetic data with known ground truth to evaluate precision and recall at different error tolerances.

**Script:** `scripts/run_precision_recall.py`
**Output:** `results/precision_recall/fig_precision_recall.png`

```bash
python scripts/run_precision_recall.py --num-reads 50000 --threads 4
```

### 3. 10x Chromium v2 Benchmark (Splitcode)
Validates Splitcode performance on 10x Chromium v2 short reads using an optimized positional extraction configuration. This was used to correct the "N/A" entries in the initial draft.

**Script:** `scripts/benchmark_splitcode_10x.py`
**Config:** `configs/splitcode/10x_v2_user.config`
**Output:** Console output (Runtime, Memory, Recovery %)

```bash
python scripts/benchmark_splitcode_10x.py
```

### 4. SPLiT-seq Concordance (Seqproc vs Split-pipe)
Calculates the Jaccard index and concordance metrics between seqproc extracted reads and the vendor's `split-pipe` output.

**Script:** `scripts/compare_splitseq.py`
**Prerequisites:** Requires `seqproc` output FASTQ and `split-pipe` output `barcode_head.fastq`.

```bash
python scripts/compare_splitseq.py \
    --seqproc /path/to/seqproc_R2.fq \
    --splitpipe /path/to/splitpipe_barcode_head.fastq
```

### 5. Sci-Seq 3 Jaccard Analysis
Performs a pairwise intersection analysis (Jaccard index) of recovered read IDs between `seqproc`, `matchbox`, and `splitcode` on the Sci-Seq 3 dataset (`SRR7827254`).

**Script:** `scripts/sciseq_jaccard_analysis.py`
**Output:** Console output (Intersection counts, Jaccard indices, Unique reads per tool)

```bash
python scripts/sciseq_jaccard_analysis.py
```

## Configurations

The analysis relies on specific configuration files for each tool:

### seqproc (`configs/seqproc/`)
- `splitseq_real.geom`: Main SPLiT-seq geometry (Anchor Relative).
- `splitseq_filter.geom`: "Locked-in" geometry with whitelist filtering (referenced in Final Report).
- `10x_longread_fwd.geom` / `_rev.geom`: Dual-pass geometry for 10x Long Read (GridION/PromethION) to handle mixed orientation.
- `sciseq3.geom`: Geometry for Sci-Seq 3 analysis.

### matchbox (`configs/matchbox/`)
- `splitseq.mb`: SPLiT-seq configuration script.
- `10x_longread.mb`: 10x Long Read configuration.

### splitcode (`configs/splitcode/`)
- `splitseq_paper.config`: Main SPLiT-seq configuration.
- `10x_v2_user.config`: Optimized positional extraction config for 10x v2.
- `sciseq3.config`: Configuration for Sci-Seq 3.

## Directory Structure

```
├── configs/                    # Tool configurations (.geom, .mb, .config)
├── data/                       # Input datasets (not included in repo)
├── docs/                       # Final reports and analysis summaries
├── paper_summaries/            # detailed daily summaries of analysis steps
├── results/
│   ├── paper_figures/          # Final Benchmark Tables & Plots
│   ├── paper_figures_revised/  # Revised SPLiT-seq specific figures
│   └── precision_recall/       # Precision-recall analysis results
└── scripts/                    # Analysis and plotting scripts
```

## Citation

If you use this analysis, please cite the seqproc paper.
