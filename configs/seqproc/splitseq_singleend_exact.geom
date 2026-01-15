# SPLiT-seq single-end with EXACT matching (no hamming tolerance)
# Fastest possible - but may miss reads with errors

umi = u[10]
bc3 = b[8]
# Exact match - no hamming distance
linker_a = anchor_relative(f[GTGGCCGATGTTTCGCATCGGCGTACGACT])
bc2 = b[8]
linker_b = anchor_relative(f[ATCCACGTGCTTGAGACTGTGG])
bc1 = b[8]
rest = r:

1{<umi><bc3><linker_a><bc2><linker_b><bc1><rest>}

-> 1{<umi><bc3><bc2><bc1>}
