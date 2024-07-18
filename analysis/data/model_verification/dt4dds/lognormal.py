import sys
import numpy as np
import dt4dds

import workflow




seqlist = dt4dds.tools.fasta_to_seqlist('./design_files.fasta')


eff_dict = {}
for seq in seqlist:
    eff_dict[seq] = 1.02-0.02*np.random.lognormal(0, 0.5)



workflow.run(eff_dict, sys.argv[1])