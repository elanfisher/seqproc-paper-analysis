primer = anchor_relative(hamming(f[CTACACGACGCTCTTCCGATCT], 3))
bc = b[16]
umi = u[12]
rest = r:
1{<primer><bc><umi><rest>}
-> 1{<bc><umi>} 2{<rest>}
