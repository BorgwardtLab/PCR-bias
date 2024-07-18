import sys
import numpy as np
from Bio import SeqIO

import workflow


eff_dict = {}
for seq in SeqIO.parse("./design_files.fasta", "fasta"):
    eff_dict[seq.id] = 1.02-0.02*np.random.lognormal(0, 0.5)



workflow.run(eff_dict, sys.argv[1])