import sys
import numpy as np
from Bio import SeqIO

import workflow


eff_dict = {}
for seq in SeqIO.parse("./design_files.fasta", "fasta"):
    eff_dict[seq.id] = np.random.uniform(0.98, 1.02)



workflow.run(eff_dict, sys.argv[1])