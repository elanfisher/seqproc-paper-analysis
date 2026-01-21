# SPLiT-seq Round 2 geometry - POSITIONAL EXTRACTION (Linkers Skipped)
# Testing user hypothesis: Extract by position, validate by list.
# Structure on R2: [NN:2][UMI:10][BC3:8][Linker1:30][BC2:8][Linker2:30][BC1:6][rest]

read1 = r:
skip2 = x[2]
umi = u[10]

# BC3: 8bp
bc3 = b[8]

# Linker1 - 30bp - SKIPPED (No validation)
skip_l1 = x[30]

# BC2: 8bp
bc2 = b[8]

# Linker2 - 30bp - SKIPPED (No validation)
skip_l2 = x[30]

# BC1: 6bp (truncated)
bc1 = b[6]
rest = r:

1{<read1>}
2{
    <skip2>
    <umi>
    # Map with mismatch using external maps passed via -a (Seq -> Seq)
    # This validates the barcode and corrects errors
    map_with_mismatch(<bc3>, $0, self, 1)
    <skip_l1>
    map_with_mismatch(<bc2>, $1, self, 1)
    <skip_l2>
    map_with_mismatch(<bc1>, $2, self, 1)
    <rest>
}

-> 1{<read1>} 2{<umi><bc3><bc2><bc1>}
