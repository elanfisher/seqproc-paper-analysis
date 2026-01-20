# 10x Chromium v2 Short Read Geometry
# R1: [Barcode:16][UMI:10]
# R2: cDNA (biological sequence)

bc = b[16]
umi = u[10]
bio = r:

# Read 1 contains barcode and UMI
1{<bc><umi>}

# Read 2 contains biological sequence
2{<bio>}

# Output: barcode, UMI, and biological read
-> 1{<bc><umi>} 2{<bio>}
