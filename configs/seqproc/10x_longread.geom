# 10x Chromium long-read barcode extraction
# Structure: [_][primer:22][barcode:16][_]
# Primer: CTACACGACGCTCTTCCGATCT (also check reverse complement)

primer = anchor_relative(hamming(f[CTACACGACGCTCTTCCGATCT], 3))
bc = b[16]
rest = r:

# Single read with primer + barcode
1{<primer><bc><rest>}

# Output: extracted barcode
-> 1{<bc>}
