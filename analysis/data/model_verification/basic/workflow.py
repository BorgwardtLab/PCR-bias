import numpy as np
import pandas as pd

PCR_EFFICIENCY = 0.9


def downsample(counts):
    # downsample counts to 1M oligos
    counts = counts / counts.sum() * 1e6
    return counts.round()



def run(eff_dict, target_folder):

    seq_ids = list(eff_dict.keys())
    n_seqs = len(seq_ids)

    # generate and save initial pool counts
    counts = downsample(np.random.lognormal(10, 0.4, n_seqs))
    counts_df = pd.DataFrame.from_dict({'seq_id': seq_ids, 'counts': counts.astype(int)})
    counts_df.to_csv(f"{target_folder}/initial_pool.csv", header=False, index=False)

    # generate and save efficiencies
    efficiencies = np.array(list(eff_dict.values())) * (1+PCR_EFFICIENCY)
    eff_df = pd.DataFrame.from_dict({'seq_id': seq_ids, 'eff': efficiencies})
    eff_df.to_csv(f"{target_folder}/eff_reference.csv", header=False, index=False)

    results_dict = {}
    # serial amplification for 6x15 cycles
    for i in range(1, 6+1):
        counts *= efficiencies**15
        counts = downsample(counts)
        results_dict[f'PCR{i}'] = counts.copy().astype(int)

    # export data
    df = pd.DataFrame.from_dict(results_dict)
    df.index = seq_ids
    df.to_csv(f"{target_folder}/abundance_by_experiment.csv", header=False)