# SPLiT-seq single-end - EXACT primer anchoring (like matchbox)
# No hamming tolerance - pure exact match

linker1 = anchor_relative(f[GTGGCCGATGTTTCGCATCGGCGTACGACT])

umi = u[10]
bc3 = b[8]
bc2 = b[8]
linker2 = x[22]
bc1 = b[8]
rest = r:

1{<umi><bc3><linker1><bc2><linker2><bc1><rest>}

-> 1{<umi><bc3><bc2><bc1>}
