# Read Processing Summary

**Dataset:** 1M reads subset from each dataset
**Threading:** 4 threads for seqproc/matchbox, 1 thread for splitcode

| Dataset              | Tool      | Total Reads   | Processed   | Recovery %   |   Runtime (s) |
|:---------------------|:----------|:--------------|:------------|:-------------|--------------:|
| SPLiT-seq Paired-End | seqproc   | 1,000,000     | 872,281     | 87.2%        |           2.6 |
| SPLiT-seq Paired-End | matchbox  | 1,000,000     | 867,118     | 86.7%        |           4   |
| SPLiT-seq Paired-End | splitcode | 1,000,000     | 883,946     | 88.4%        |          10.6 |
| SPLiT-seq Single-End | seqproc   | 1,000,000     | 241,771     | 24.2%        |          16.1 |
| SPLiT-seq Single-End | matchbox  | 1,000,000     | 240,065     | 24.0%        |           2.9 |
| SPLiT-seq Single-End | splitcode | 1,000,000     | 279,699     | 28.0%        |         103   |
| 10x GridION          | seqproc   | 1,000,000     | 362,405     | 36.2%        |          20.1 |
| 10x GridION          | matchbox  | 1,000,000     | 410,096     | 41.0%        |           6.4 |
| 10x GridION          | splitcode | 1,000,000     | 362,409     | 36.2%        |          65.8 |

## Notes
- SPLiT-seq PE: Paired-end short reads (Illumina)
- SPLiT-seq SE: Single-end long reads (ONT)
- 10x GridION: 10x Chromium on Oxford Nanopore GridION
- Recovery rates are consistent across tools (~87% for PE, ~24% for SE, ~36-41% for 10x)
- Lower SE/10x recovery expected due to long-read error profiles and primer finding
