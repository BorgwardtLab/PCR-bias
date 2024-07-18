import dt4dds
import pandas as pd
import shutil
import gc

PCR_EFFICIENCY = 0.9

def run(eff_dict, target_folder):

    dt4dds.config.show_progressbars = False
    dt4dds.config.enable_multiprocessing = False
    dt4dds.config.n_processes = 1
    dt4dds.default_logging()


    primers_0 = ["ACACGACGCTCTTCCGATCT", "AGACGTGTGCTCTTCCGATCT"]
    primers_2 = ["AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT", "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"]

    # assign efficiency properties
    dt4dds.properties.AmplificationEfficiency._efficiency_cache = eff_dict
    def disabled(*args):
        raise NotImplementedError('Amplification efficiency distribution disabled')
    dt4dds.properties.AmplificationEfficiency._efficiency_distribution = disabled
    df = pd.DataFrame.from_dict(eff_dict, columns=["eff"], orient="index")
    df.eff = df.eff * (1+PCR_EFFICIENCY)
    df.to_csv(f"{target_folder}/eff_reference.csv", header=False)

    # common settings for PCR
    pcr_settings = dt4dds.settings.defaults.PCR(
        primers=primers_0,
        template_volume=5,
        volume=20,
        efficiency_mean=PCR_EFFICIENCY,
        polymerase_basesubstitutionrate=0,
        n_cycles=15
    )


    # common sequencing workflow for bunnies
    def sequence(seqpool, name):

        # samples dilution
        pool = seqpool.sample_by_volume(5)

        # sequencing PCR with indexed primers
        pcr = dt4dds.processes.PCR(pcr_settings(primers=primers_2))
        pool = pcr.process(pool)

        # dilution to 50 pM * 20 uL
        pool = pool.sample_by_moles(0.05*20e-6)

        # sequencing
        sbs_sequencing = dt4dds.processes.SBSSequencing(
            dt4dds.settings.defaults.iSeq100(
                output_directory=f"{target_folder}/{name}/",
                n_reads=1000000,
                substitution_rates=[0, 0]
            )
        )
        sbs_sequencing.process(pool)
        shutil.copyfile('./design_files.fasta', f"{target_folder}/{name}/design_files.fasta")





    # read sequences from reference
    seq_list = dt4dds.tools.fasta_to_seqlist('./design_files.fasta')

    # use synthesis settings as passed
    array_synthesis = dt4dds.processes.ArraySynthesis(dt4dds.settings.defaults.ArraySynthesis_Twist(
        oligo_distribution_type='lognormal',
        oligo_distribution_params={'mean': 0, 'sigma': 0.30},
        deletion_rate=0,
        insertion_rate=0,
        substition_rate=0,
    ))
    array_synthesis.process(seq_list)

    # perform synthesis
    pool = array_synthesis.sample_by_mass(1.0048*1/50)

    # save initial distribution
    pool.save_as_csv(f"{target_folder}/initial_pool.csv")

    # attach primers
    pool = dt4dds.generators.attach_primers_to_pool(pool, *primers_0)
    pool = pool.sample_by_mass(1.0048*1/50/10)
    pool.volume = 20


    # 
    # PCR0 pool
    # 
    sequence(pool, f'PCR1')
    print("PCR1 done")
    gc.collect()


    # 
    # multiple amplifications experiments
    # 
    for i in range(1, 6):

        # amplification PCR
        pcr = dt4dds.processes.PCR(pcr_settings)
        pool = pcr.process(pool)
        
        # dilute
        pool = pool.sample_by_volume(1)
        pool.volume = 100
        pool = pool.sample_by_volume(1)
        pool.volume = 38

        
        # sequencing
        sequence(pool, f'PCR{i+1}')
        print(f'PCR{i} done')
        gc.collect()