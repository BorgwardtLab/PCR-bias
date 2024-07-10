from captum.attr import DeepLift
import torch
import pandas as pd
import numpy as np
from utils.utils import representation
import matplotlib.pyplot as plt
import logomaker
import warnings
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import average_precision_score, roc_auc_score
from model import CNN
import pickle
from pathlib import Path
from utils.motif_utils import (
    find_most_significant_kmer_in_sequence,
    one_hot_to_DNA,
    pairwise_distance,
    tsne_kmeans_search_silhouette,
    create_pwm_for_cluster,
    convolve_pwm_with_sequences,
    perform_chi_squared_tests
)

warnings.filterwarnings("ignore")

class CLAMAnalysis:
    def __init__(self, filename: str, threshold: int):
        self.filename = filename
        self.threshold = threshold
        self.motif_save_path = f"Motifs/{filename}/{threshold}/"
        self.model_path = f"CNN/model/internal/{filename}/{threshold}/"
        self.attr_save_path = f"Feature_attribution/{filename}/{threshold}/"
        Path(self.attr_save_path).mkdir(parents=True, exist_ok=True)
        Path(self.motif_save_path).mkdir(parents=True, exist_ok=True)

    def get_data_with_q(self, source, threshold):
        seqs = pd.read_pickle(f"data/{source}/bad_seqs_{threshold}.pkl")
        rest_seqs = seqs["rest"]["sequence"]
        bot_seqs = seqs["bottom"]["sequence"]
        q = pd.concat([seqs["bottom"].eff, seqs["rest"].eff])
        seqs = np.hstack([bot_seqs, rest_seqs])
        labels = [1] * len(bot_seqs) + [0] * len(rest_seqs)

        X_src = representation(seqs, with_reverse=False)
        X_src = torch.tensor(X_src)
        y_src = torch.tensor(labels)
        return X_src, y_src, q

    def plot_pwm_logo(self, pwm, cluster_idx):
        plt.figure(figsize=(12, 6))
        nucleotides = ["A", "C", "G", "T"]
        pwm_df = pd.DataFrame(pwm.T, columns=nucleotides)

        pwm_df += 1e-9
        entropy = -np.sum(pwm_df.apply(lambda x: x * np.log2(x)), axis=1)
        information_content = 2 - entropy
        scaled_pwm = pwm_df.mul(information_content, axis=0)
        logo = logomaker.Logo(scaled_pwm)
        plt.ylim([0, 2])
        plt.xlabel("Position")
        plt.ylabel("Information Content (bits)")
        plt.title(f"Motif Cluster {cluster_idx} - PWM Logo")
        plt.tight_layout()
        plt.savefig(
            f"{self.motif_save_path}/motif_cluster_length{pwm.shape[1]}_cluster{cluster_idx}.png",
            dpi=200,
        )
        plt.close()

    def feature_attribution(self):
        X_src, y_src, q_src = self.get_data_with_q(self.filename, self.threshold)
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=1)
        all_attributions, all_qs, all_labels, all_seqs = [], [], [], []
        # Load the model and perform feature attribution
        for fold_count, (train_idx, test_idx) in enumerate(cv.split(X_src, y_src)):
            X_train, X_test = X_src[train_idx], X_src[test_idx]
            y_train, y_test = y_src[train_idx], y_src[test_idx]
            q_train, q_test = q_src[train_idx], q_src[test_idx]

            model = CNN.load_from_checkpoint(f"{self.model_path}/best_model_fold_{fold_count}.ckpt")
            dl = DeepLift(model)  # can be replaced by any kind of feature attribution method
            attributions = dl.attribute(X_test.float(), target=1)
            all_qs.append(q_test)
            all_attributions.append(attributions.detach().numpy())
            all_labels.append(y_test)
            all_seqs.append(X_test)

        all_seqs = np.vstack(all_seqs)
        all_qs = pd.concat(all_qs)
        all_attributions = np.vstack(all_attributions)
        all_labels = np.hstack(all_labels)

        with open(f"{self.attr_save_path}/attributions.pkl", "wb") as fp:
            pickle.dump([all_attributions, all_seqs, all_labels, all_qs], fp)

    def replace_motifs(self, X_tar, significant_motifs, top_k):
        sequence_length = X_tar.shape[1]
        for k, (p_values, pwm) in enumerate(significant_motifs):
            pwm = torch.tensor(pwm)
            pwm_length = pwm.shape[1]
            avg_subseq = torch.full((pwm_length, 4), 0.25)
            for j, seq in enumerate(X_tar):
                for i in range(sequence_length - pwm_length + 1):
                    window = X_tar[j, i:i + pwm_length, :]
                    score = torch.sum(pwm * window.T) / pwm_length
                    if score >= 0.5:
                        X_tar[j, i: i + pwm_length, :] = avg_subseq
            if k == top_k:
                break
        return X_tar

    def motif_substituition_internal(self, target, threshold):
        significant_motifs = pd.read_pickle(f"Motifs/{target}/{threshold}/significantly_enriched_motifs.pkl").sort(
            key=lambda x: x[2])
        X_src, y_src, q_src =  self.get_data_with_q(target, threshold)
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=1)
        metrics_summary = pd.DataFrame(columns=["test_auprc", "test_auroc"],
                                       index=pd.MultiIndex.from_product([range(5), range(len(significant_motifs))]))

        for fold_count, (train_idx, test_idx) in enumerate(cv.split(X_src, y_src)):
            X_train, X_test = X_src[train_idx], X_src[test_idx]
            y_train, y_test = y_src[train_idx], y_src[test_idx]

            model = CNN.load_from_checkpoint(f"CNN/model/internal/{target}/{threshold}/best_model_fold_{fold_count}.ckpt")
            for top_k in range(len(significant_motifs)):
                X_test_sub = self.replace_motifs(X_test.clone(), significant_motifs, top_k)
                pred_test= model(X_test_sub.float())[:, 1]
                auprc = average_precision_score(y_test.numpy(), pred_test.detach().numpy())
                auroc = roc_auc_score(y_test.detach().numpy(), pred_test.detach().numpy())
                metrics_summary.loc[fold_count, top_k] = [auprc, auroc]
        metrics_summary.to_csv(f"CNN/results/internal/{target}/{threshold}/motif_substitution_internal.csv")

    def motif_substituition_external(self, source, target, threshold):
        significant_motifs = pd.read_pickle(f"Motifs/{target}/{threshold}/significantly_enriched_motifs.pkl").sort(
            key=lambda x: x[2])
        X_test, y_test, q_test = self.get_data_with_q(target, threshold)
        model = CNN.load_from_checkpoint(f"CNN/model/external/{source}/{threshold}/best_model.ckpt")
        metrics_summary = pd.DataFrame(columns=["test_auprc", "test_auroc"], index=range(len(significant_motifs)))

        for top_k in range(len(significant_motifs)):
            X_test_sub = self.replace_motifs(X_test.clone(), significant_motifs, top_k)
            pred_test = model(X_test_sub.float())[:, 1]
            auprc = average_precision_score(y_test.numpy(), pred_test.detach().numpy())
            auroc = roc_auc_score(y_test.detach().numpy(), pred_test.detach().numpy())
            metrics_summary.loc[top_k] = [auprc, auroc]
        metrics_summary.to_csv(f"CNN/results/internal/{source}/{threshold}/motif_substitution_{target}.csv")

    def PWM_construction(self):
        seqs_attr, seqs_onehot, seqs_label, seqs_q = np.load(f"{self.attr_save_path}/attributions.pkl", allow_pickle=True)
        seqs_attr = seqs_attr[np.argsort(seqs_q.tolist())]
        seqs_label = seqs_label[np.argsort(seqs_q.tolist())]
        bot_seqs_attr = seqs_attr[np.where(seqs_label == 1)[0]]

        window_sizes = np.arange(4, 13)
        pwms_per_window_size = {}
        # clustering on all windows sizes
        for window_size_idx, window_size in enumerate(window_sizes):
            motifs_per_window_size = []
            for seq_idx in range(bot_seqs_attr.shape[0]):
                kmer_start, _ = find_most_significant_kmer_in_sequence(
                    bot_seqs_attr[seq_idx], window_size=window_size, stride=1
                )
                if kmer_start is not None:
                    kmer = bot_seqs_attr[seq_idx, kmer_start: kmer_start + window_size, :]
                    alphabetic_kmer = one_hot_to_DNA((kmer != 0).astype(int))
                    if len(alphabetic_kmer) == window_size:
                        motifs_per_window_size.append(alphabetic_kmer)

            motifs_per_window_size = pd.Series(motifs_per_window_size)
            unique_motifs = pd.Series(motifs_per_window_size.unique())
            distance_matrix = pd.DataFrame(
                pairwise_distance(unique_motifs),
                index=unique_motifs,
                columns=unique_motifs,
            )
            sample_weight = motifs_per_window_size.value_counts().loc[distance_matrix.index]
            clusters, embeddings_2d = tsne_kmeans_search_silhouette(
                distance_matrix, sample_weight
            )
            embeddings_2d = pd.DataFrame(
                embeddings_2d, columns=["dim 1", "dim 2"], index=distance_matrix.index
            )
            embeddings_2d["occurrence"] = motifs_per_window_size.value_counts().loc[
                embeddings_2d.index
            ]
            embeddings_2d["cluster"] = clusters

            embeddings_2d.to_csv(f"{self.motif_save_path}/motif_cluster_length{window_size}.csv")
            pwm = embeddings_2d[["cluster", "occurrence"]]
            pwm["sequence"] = pwm.index
            pwms = {}
            for cluster in pwm["cluster"].unique():
                cluster_df = pwm[pwm["cluster"] == cluster]
                pwms[cluster] = create_pwm_for_cluster(cluster_df)
            pwms_per_window_size[window_size] = pwms

        with open(f"{self.motif_save_path}/pwms_by_window_size_whole_sequence.pkl", "wb") as fp:
            pickle.dump(pwms_per_window_size, fp)
        return seqs_onehot, seqs_label, pwms_per_window_size, window_sizes

    def statistical_test(self):
        seqs_onehot, seqs_label, pwms_per_window_size, window_sizes = self.PWM_construction()

        seqs_positive = seqs_onehot[np.array(seqs_label) == 1]
        seqs_negative = seqs_onehot[np.array(seqs_label) == 0]
        p_values_per_pwm = []

        for window_size_idx, window_size in enumerate(window_sizes):
            for cluster_idx, pwm in pwms_per_window_size[window_size].items():
                counts_positive = convolve_pwm_with_sequences(pwm, seqs_positive)
                counts_negative = convolve_pwm_with_sequences(pwm, seqs_negative)
                try:
                    chi2, p_value = perform_chi_squared_tests(counts_positive, counts_negative)
                    p_values_per_pwm.append([pwm, chi2, p_value, cluster_idx, counts_positive, counts_negative])
                except ValueError:
                    pass

        corrected_alpha = 0.05 / len(p_values_per_pwm)
        p_values_per_pwm = [_ for _ in p_values_per_pwm if _[2] < corrected_alpha]
        with open(f"{self.motif_save_path}/significantly_enriched_motifs.pkl", "wb") as fp:
            pickle.dump(p_values_per_pwm, fp)

        for pwm, _, p_value, cluster_idx, _, _ in p_values_per_pwm:
            self.plot_pwm_logo(pwm, cluster_idx)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--filename", type=str, required=True)
    parser.add_argument("--threshold", type=str, required=True)
    args = parser.parse_args()

    clam_analysis = CLAMAnalysis(args.filename, args.threshold)
    clam_analysis.feature_attribution()
    clam_analysis.statistical_test()
