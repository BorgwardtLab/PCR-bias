import numpy as np
import pandas as pd
import scipy.stats
import scipy.sparse
import scipy.sparse.linalg


def read_csvs(dir, keep_all_sequences=False):
    df = pd.read_csv(f"{dir}/abundance_by_experiment.csv", dtype={'seq_id': str}, index_col="seq_id")
    
    for col in df.columns:
        print(f"{col}: {df[col].sum()} total reads, {len(df[col])-np.count_nonzero(df[col])} missing sequences")

    # remove sequences with reads in only one or no experiment
    seqs_by_occurrence = df[df > 0].count(axis=1)
    print(f"{seqs_by_occurrence.value_counts().get(1,0)} sequences with reads in only one experiment, {seqs_by_occurrence.value_counts().get(0,0)} sequences with no reads at all.")
    if not keep_all_sequences: df = df[seqs_by_occurrence > 1]

    # normalize each experiment by the mean number of reads
    df = df.div(df.mean(axis=0), axis=1)

    return df


def fit_parameters(observations_df, experimental_parameters):
    # melt the dataframe
    n_seqs = len(observations_df)
    wide_df = observations_df.sort_index()
    wide_df['seq_id'] = wide_df.index
    df = pd.melt(wide_df, "seq_id", value_vars=wide_df.columns, var_name="exp", value_name="x")
    df = df.sort_values(by=['seq_id', 'exp']).reset_index(drop=True)

    # add n_cycles column
    df['n_cycles'] = 0
    for exp in df.exp.unique():
        df.loc[df.exp == exp, 'n_cycles'] = experimental_parameters[exp]

    # prepare the least-squares matrix with all the data
    cycles = np.array([value for key, value in sorted(experimental_parameters.items())])
    initials = np.ones((len(experimental_parameters),))
    A = scipy.sparse.lil_matrix((len(cycles)*n_seqs, 2*n_seqs), dtype=float)
    for i in range(n_seqs):
        A[len(cycles)*i:len(cycles)*(i+1), i] = cycles
        A[len(cycles)*i:len(cycles)*(i+1), n_seqs+i] = initials

    # mask the least-squares matrix and the output vector to only include nonzero values
    nonzero_mask = (df.x != 0).values
    A = A.tocsr()[nonzero_mask]
    b = np.log(df.x[nonzero_mask].values)

    # perform sparse least squares for log-model
    rough_result = scipy.sparse.linalg.lsqr(
        A, 
        b, 
        show=True, 
        atol=1e-12, 
        btol=1e-12, 
        iter_lim=1000
    )

    # from the results, extract non-log parameters of the model
    x0 = np.exp(rough_result[0][n_seqs:])
    x0 /= x0.mean()
    eff = np.exp(rough_result[0][:n_seqs])
    eff /= eff.mean()

    # prepare dataframe with parameter results
    return pd.DataFrame.from_dict({'seq_id': wide_df.seq_id, 'x0': x0, 'eff': eff}).set_index('seq_id')


def assess_fit(observations_df, parameter_df, experimental_parameters):
    
    # create a dataframe with model predictions
    cols = {}
    for col, n_cyc in experimental_parameters.items():
        cols[col] = parameter_df.x0 * (parameter_df.eff**n_cyc)
    model_df = pd.DataFrame(cols)
    model_df = model_df.div(model_df.mean(axis=0), axis=1)
    model_df = pd.melt(model_df.reset_index(), id_vars="seq_id", var_name="exp", value_name="model")
    
    # create a dataframe with a simple baseline model
    cols = {}
    for col, n_cyc in experimental_parameters.items():
        cols[col] = observations_df.mean(axis=1)
    baseline_df = pd.DataFrame(cols)
    baseline_df = pd.melt(baseline_df.reset_index(), id_vars="seq_id", var_name="exp", value_name="baseline")
    
    # melt the observations dataframe
    data_df = pd.melt(observations_df.reset_index(), id_vars="seq_id", var_name="exp", value_name="true")

    # merge both models with the "true" observations
    all_df = pd.merge(data_df, model_df, on=["seq_id", "exp"])
    all_df = pd.merge(all_df, baseline_df, on=["seq_id", "exp"])
    
    # compile fit assessment
    d = {}

    # calculate RSS and explained variance
    d["rss_tot"] = sum((all_df.baseline.mean()-all_df.true)**2)
    d["rss_baseline"] = sum((all_df.baseline-all_df.true)**2)
    d["rss_model"] = sum((all_df.model-all_df.true)**2)
    d["expvar_baseline"] = 1 - d["rss_baseline"] / d["rss_tot"]
    d["expvar_model"] = 1 - d["rss_model"] / d["rss_tot"]

    # calculate F-statistic
    d["f_df1"] = 2*len(parameter_df) - len(parameter_df)
    d["f_df2"] = len(data_df) - 2*len(parameter_df)
    if d["f_df2"] <= 0:
        d["f_value"] = None
        d["f_pvalue"] = None
    else:
        d["f_value"] = ((d["rss_baseline"]-d["rss_model"])/d["f_df1"]) / (d["rss_model"]/d["f_df2"])
        d["f_pvalue"] = 1 - scipy.stats.f.cdf(d["f_value"], dfn=d["f_df1"], dfd=d["f_df2"])

    # return merged df and assessment
    return all_df, d