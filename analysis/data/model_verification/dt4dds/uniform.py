import sys
import numpy as np
import dt4dds

import workflow




seqlist = dt4dds.tools.fasta_to_seqlist('./design_files.fasta')


eff_dict = {}
for seq in seqlist:
    eff_dict[seq] = np.random.uniform(0.98, 1.02)



workflow.run(eff_dict, sys.argv[1])