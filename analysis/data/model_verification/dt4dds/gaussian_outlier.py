import sys
import numpy as np
import dt4dds

import workflow




seqlist = dt4dds.tools.fasta_to_seqlist('./design_files.fasta')


eff_dict = {}
for seq in seqlist:
    if np.random.random() < 0.02:
        eff_dict[seq] = np.random.uniform(0.9, 0.98)
    else:
        eff_dict[seq] = np.random.normal(1.0, 0.0083)



workflow.run(eff_dict, sys.argv[1])