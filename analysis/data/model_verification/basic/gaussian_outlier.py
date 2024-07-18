import sys
import numpy as np
from Bio import SeqIO

import workflow


eff_dict = {}
for seq in SeqIO.parse("./design_files.fasta", "fasta"):
    if np.random.random() < 0.02:
        eff_dict[seq.id] = np.random.uniform(0.9, 0.98)
    else:
        eff_dict[seq.id] = np.random.normal(1.0, 0.0083)



workflow.run(eff_dict, sys.argv[1])