primer_rc = anchor_relative(hamming(f[AGATCGGAAGAGCGTCGTGTAG], 3))
bc_rc = b[16]
rest = r:

1{<rest><bc_rc><primer_rc>}

-> 1{revcomp(<bc_rc>)}
