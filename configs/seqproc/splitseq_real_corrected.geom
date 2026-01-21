# SPLiT-seq Round 2 geometry for REAL data (SRR6750041)
# With Whitelist Correction (map_with_mismatch)
# Structure on R2: [NN:2][UMI:10][BC3:8][Linker1:30][BC2:8][Linker2:30][BC1:6][rest]

read1 = r:
skip2 = x[2]
umi = u[10]

# BC3: 8bp, whitelist_v1.txt, 1 mismatch
bc3 = b[8]

# Linker1 - 30bp - Use relative anchor with moderate strictness (Hamming 5)
l1 = anchor_relative(hamming(f[GTGGCCGCTGTTTCGCATCGGCGTACGACT], 5))

# BC2: 8bp, whitelist_v1.txt, 1 mismatch
bc2 = b[8]

# Linker2 - 30bp - Use relative anchor with moderate strictness (Hamming 5)
l2 = anchor_relative(hamming(f[ATCCACGTGCTTGAGAGGCCAGAGCATTCG], 5))

# BC1: 6bp (truncated), whitelist_v1_6bp.txt, 1 mismatch
bc1 = b[6]
rest = r:

1{<read1>}
2{
    <skip2>
    <umi>
    map_with_mismatch(<bc3>, $0, self, 1)
    <l1>
    map_with_mismatch(<bc2>, $1, self, 1)
    <l2>
    map_with_mismatch(<bc1>, $2, self, 1)
    <rest>
}

-> 1{<read1>} 2{<umi><bc3><bc2><bc1>}
