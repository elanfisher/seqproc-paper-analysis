# SPLiT-seq single-end with PRIMER ANCHORING
# Strategy: Find Linker1 once, then use fixed offsets for all barcodes
#
# Structure: [UMI:10][BC3:8][Linker1:30][BC2:8][Linker2:22][BC1:8][rest]
#
# Anchor on Linker1, extract everything relative to it

# First, anchor on Linker1 (the main anchor point)
linker1 = anchor_relative(hamming(f[GTGGCCGATGTTTCGCATCGGCGTACGACT], 5))

# Before linker1: UMI (10bp) + BC3 (8bp) = 18bp before linker1
umi = u[10]
bc3 = b[8]

# After linker1: BC2 (8bp) + Linker2 (22bp) + BC1 (8bp)
bc2 = b[8]
linker2 = x[22]  # Skip 22bp (don't search, fixed offset)
bc1 = b[8]
rest = r:

# Parse structure: anchor finds linker1, then fixed positions before/after
1{<umi><bc3><linker1><bc2><linker2><bc1><rest>}

# Output: extracted barcodes
-> 1{<umi><bc3><bc2><bc1>}
