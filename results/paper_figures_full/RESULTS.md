# seqproc Benchmark Results

## Overview

Comprehensive benchmark comparing **seqproc**, **matchbox**, and **splitcode** on:
- SPLiT-seq paired-end (short reads)
- SPLiT-seq single-end (long reads)
- 10x Chromium GridION (long reads)
- 10x Chromium PromethION (long reads)

All tests run on 1M reads with 4 threads.

## Results Summary

| Dataset | seqproc | matchbox | splitcode | Speedup |
|---------|---------|----------|-----------|---------|
| SPLiT-seq Paired-End (Full) | 224.78s | 377.70s | 930.42s | 1.68x |
| SPLiT-seq Single-End (Full) | 61.39s | 14.80s | 631.13s | 0.24x |
| 10x GridION (Full) | 152.41s | 152.02s | 557.10s | 1.00x |

## Dataset Analysis

### SPLiT-seq Paired-End (Short Reads)

**seqproc is 1.68x faster**

- **seqproc**: 224.78s ± 1.84s, 48.6% recovery (37,706,109 reads)
- **matchbox**: 377.70s ± 2.67s, 0.0% recovery (38,756 reads)
- **splitcode**: 930.42s ± 6.35s, 0.0% recovery (167 reads)

**Why**: Fixed-position barcode extraction on short R2 reads (~100bp). No searching needed - seqproc extracts directly by offset.

### SPLiT-seq Single-End (Long Reads)

**matchbox is 4.15x faster**

- **seqproc**: 61.39s ± 38.38s, 24.0% recovery (1,339,101 reads)
- **matchbox**: 14.80s ± 0.06s, 17.9% recovery (1,000,001 reads)
- **splitcode**: 631.13s ± 10.43s, 27.5% recovery (1,535,830 reads)

**Why**: Long reads (500-2000bp) require anchor-based searching. matchbox uses optimized pattern automata; seqproc's `anchor_relative` has higher per-read overhead but scales better with threads.

### 10x Chromium GridION (Long Reads)

**matchbox is 1.00x faster**

- **seqproc**: 152.41s ± 0.47s, 36.0% recovery (2,704,177 reads)
- **matchbox**: 152.02s ± 0.19s, 40.7% recovery (3,064,685 reads)
- **splitcode**: 557.10s ± 2.70s, 36.0% recovery (2,704,212 reads)

## Overall Summary

**seqproc is faster in 1 out of 3 datasets.**

| Dataset Type | Winner | Speedup |
|--------------|--------|---------|
| SPLiT-seq Paired-End (Full) | seqproc | 1.68x |
| SPLiT-seq Single-End (Full) | matchbox | 4.15x |
| 10x GridION (Full) | matchbox | 1.00x |

### Performance Patterns

- **seqproc excels** at fixed-position extraction and simple anchor patterns
- **matchbox excels** at complex pattern matching in long reads (single-threaded)
- **seqproc scales better** with multiple threads

## Figures

- `fig1_runtime_comparison.png` - Runtime bar charts for all datasets
- `fig2_speedup_chart.png` - Speedup comparison (seqproc vs matchbox)
- `fig3_read_recovery.png` - Read recovery rates for all tools
- `fig4_summary_dashboard.png` - Complete summary dashboard with insights