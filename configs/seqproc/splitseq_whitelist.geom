# SPLiT-seq with WHITELIST-BASED barcode matching (like splitcode)
# Uses search_whitelist to find and REPLACE barcodes with canonical sequences
#
# Structure on R2: [NN:2][UMI:10][BC3:8][Linker1:30][BC2:8][Linker2:16][BC1:8][rest]
#
# Whitelist files contain one barcode per line
# search_whitelist(b[n], file, hamming_dist) searches for barcode from file
# and replaces with the matching canonical barcode

read1 = r:
skip2 = x[2]
umi = u[10]

# BC3 with whitelist matching (2 mismatches allowed)
bc3 = search_whitelist(b[8], "configs/seqproc/splitseq_bc_whitelist.txt", 2)

# Linker1 - 30bp anchor
l1 = anchor_relative(hamming(f[GTGGCCGCTGTTTCGCATCGGCGTACGACT], 9))

# BC2 with whitelist matching
bc2 = search_whitelist(b[8], "configs/seqproc/splitseq_bc_whitelist.txt", 2)

# Linker2 - 16bp anchor  
l2 = anchor_relative(hamming(f[ATCCACGTGCTTGAGA], 5))

# BC1 with whitelist matching
bc1 = search_whitelist(b[8], "configs/seqproc/splitseq_bc_whitelist.txt", 2)

rest = r:

1{<read1>}
2{<skip2><umi><bc3><l1><bc2><l2><bc1><rest>}

-> 1{<read1>} 2{<umi><bc3><bc2><bc1>}
