# SPLiT-seq single-end OPTIMIZED - reduced hamming tolerance
# Trade accuracy for speed

umi = u[10]
bc3 = b[8]
# Reduced from hamming 9 to hamming 4 (13% vs 30% tolerance)
linker_a = anchor_relative(hamming(f[GTGGCCGATGTTTCGCATCGGCGTACGACT], 4))
bc2 = b[8]
# Reduced from hamming 7 to hamming 3 (14% vs 32% tolerance)
linker_b = anchor_relative(hamming(f[ATCCACGTGCTTGAGACTGTGG], 3))
bc1 = b[8]
rest = r:

1{<umi><bc3><linker_a><bc2><linker_b><bc1><rest>}

-> 1{<umi><bc3><bc2><bc1>}
