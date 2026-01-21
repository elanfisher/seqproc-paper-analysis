# SPLiT-seq with FAST barcode replacement using pre-computed mismatch variants
# Uses exact hash lookup instead of hamming distance matching
# For REAL data (SRR6750041) - GCT variant linkers
# Structure on R2: [NN:2][UMI:10][BC3:8][Linker1:30][BC2:8][Linker2:16][BC1:8][rest]

read1 = r:
skip2 = x[2]
umi = u[10]

# BC3 with exact lookup (mismatch variants pre-computed in expanded TSV)
bc3_def = b[8]

# Linker1 - 30bp with hamming distance 6 (~20% error rate)
l1 = anchor_relative(hamming(f[GTGGCCGCTGTTTCGCATCGGCGTACGACT], 6))

# BC2 with exact lookup
bc2_def = b[8]

# Linker2 - 16bp with hamming distance 3 (~20% error rate)
l2 = anchor_relative(hamming(f[ATCCACGTGCTTGAGA], 3))

# BC1 with exact lookup
bc1_def = b[8]
rest = r:

1{<read1>}
2{
    <skip2>
    <umi>
    map(<bc3_def>, $0, self)
    <l1>
    map(<bc2_def>, $1, self)
    <l2>
    map(<bc1_def>, $2, self)
    <rest>
}

-> 1{<read1>} 2{<umi><bc3_def><bc2_def><bc1_def>}
