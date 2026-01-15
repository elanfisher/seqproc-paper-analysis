# SPLiT-seq geometry for SINGLE-END long-read data (SRR13948564)
# Structure: [UMI:10][BC3:8][Linker1:30][BC2:8][Linker2:22][BC1:8][rest]
#
# Single-read mode: file1 only, file2 uses same file automatically
# Uses anchor_relative to find linkers with fuzzy matching (hamming distance)

read1 = r:
umi = u[10]
bc3 = b[8]
linker_a = anchor_relative(hamming(f[GTGGCCGATGTTTCGCATCGGCGTACGACT], 9))
bc2 = b[8]
linker_b = anchor_relative(hamming(f[ATCCACGTGCTTGAGACTGTGG], 7))
bc1 = b[8]
rest = r:

# Read 1 is dummy (same file), Read 2 has the barcode structure
1{<read1>}
2{<umi><bc3><linker_a><bc2><linker_b><bc1><rest>}

# Output: R1 unchanged, R2 with extracted UMI+barcodes
-> 1{<read1>} 2{<umi><bc3><bc2><bc1>}
