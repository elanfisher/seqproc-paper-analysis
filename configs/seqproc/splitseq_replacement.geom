# SPLiT-seq with barcode REPLACEMENT using map_with_mismatch
# For REAL data (SRR6750041) - GCT variant linkers
# CORRECTED Structure on R2: [NN:2][UMI:8][BC3:8][Linker1:30][BC2:8][Linker2:30][BC1:8][rest]
# Note: UMI is 8bp (not 10bp), Linker2 is 30bp (not 16bp)

read1 = r:
skip2 = x[2]
umi = u[8]

# BC3 with whitelist replacement (2 mismatch tolerance)
bc3_def = b[8]

# Linker1 - 30bp with hamming distance 2 (same as splitcode)
l1 = anchor_relative(hamming(f[GTGGCCGCTGTTTCGCATCGGCGTACGACT], 2))

# BC2 with whitelist replacement
bc2_def = b[8]

# Linker2 - 30bp with hamming distance 2 (CORRECTED - was incorrectly 16bp)
l2 = anchor_relative(hamming(f[ATCCACGTGCTTGAGAGGCCAGAGCATTCG], 2))

# BC1 with whitelist replacement  
bc1_def = b[8]
rest = r:

1{<read1>}
2{
    <skip2>
    <umi>
    map_with_mismatch(<bc3_def>, $0, self, 2)
    <l1>
    map_with_mismatch(<bc2_def>, $1, self, 2)
    <l2>
    map_with_mismatch(<bc1_def>, $2, self, 2)
    <rest>
}

-> 1{<read1>} 2{<umi><bc3_def><bc2_def><bc1_def>}
