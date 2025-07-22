from captum.attr import DeepLift
import torch
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logomaker
import pickle
from pathlib import Path

from utils.motif_utils import (
    find_most_significant_kmer_in_sequence,
    one_hot_to_DNA,
    pairwise_distance,
    tsne_kmeans_search_silhouette,
    create_pwm_for_cluster,
    convolve_pwm_with_sequences,
    perform_chi_squared_tests,
    representation,
)
import os
os.environ["PYTHONWARNINGS"] = "ignore"
import warnings
warnings.filterwarnings("ignore")


class CluMo:
    def __init__(self, filename: str, threshold: int):
        self.filename = filename
        self.threshold = threshold
        self.motif_save_path = f"CNN/motifs/{filename}/{threshold}/"
        self.model_path = f"CNN/models/"
        Path(self.motif_save_path).mkdir(parents=True, exist_ok=True)

    def get_data_with_q(self, source, threshold):
        seqs = pd.read_pickle(f"Data/{source}/bad_seqs_{threshold}.pkl")
        rest_seqs = seqs["rest"]["sequence"]
        bot_seqs = seqs["bottom"]["sequence"]
        q = pd.concat([seqs["bottom"].eff, seqs["rest"].eff])
        seqs = np.hstack([bot_seqs, rest_seqs])
        labels = [1] * len(bot_seqs) + [0] * len(rest_seqs)

        X_src = representation(seqs, with_reverse=False)
        X_src = torch.tensor(X_src)
        y_src = torch.tensor(labels)
        return X_src, y_src, q

    def plot_pwm_logo(self, pwm, p_value, cluster_idx):
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
        plt.title(f"PWM Logo, p-value: {p_value:.2e}")
        plt.tight_layout()
        plt.savefig(
            f"{self.motif_save_path}/motif_cluster_length{pwm.shape[1]}_No{cluster_idx}.png",
            dpi=200,
        )

        plt.close()

    def feature_attribution(self):
        seqs_onehot, seqs_label, seqs_q = self.get_data_with_q(
            self.filename, self.threshold
        )

        seqs_attr = []
        model = torch.load(
            "CNN/models/{}_Filtered_{}_best_model.pt".format(
                self.filename, self.threshold
            ),
            map_location="cpu",
        )
        dl = DeepLift(model)
        attributions = dl.attribute(seqs_onehot.float(), target=1)
        seqs_attr.append(attributions.detach().numpy())
        seqs_attr = np.vstack(seqs_attr)
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
                    kmer = bot_seqs_attr[
                        seq_idx, kmer_start : kmer_start + window_size, :
                    ]
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
            sample_weight = motifs_per_window_size.value_counts().loc[
                distance_matrix.index
            ]
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
            pwm = embeddings_2d[["cluster", "occurrence"]]
            pwm["sequence"] = pwm.index
            pwms = {}
            for cluster in pwm["cluster"].unique():
                cluster_df = pwm[pwm["cluster"] == cluster]
                pwms[cluster] = create_pwm_for_cluster(cluster_df)
            pwms_per_window_size[window_size] = pwms
        return seqs_onehot, seqs_label, pwms_per_window_size, window_sizes

    def motif_plot(self, only_visualization=True):
        if only_visualization:
            p_values_per_pwm = pd.read_pickle(
                f"{self.motif_save_path}/significantly_enriched_motifs.pkl"
            )
            corrected_alpha = 0.05 / len(p_values_per_pwm)
            p_values_per_pwm = [_ for _ in p_values_per_pwm if _[2] < corrected_alpha]
        else:
            (
                seqs_onehot,
                seqs_label,
                pwms_per_window_size,
                window_sizes,
            ) = self.feature_attribution()

            seqs_positive = seqs_onehot[np.array(seqs_label) == 1].numpy()
            seqs_negative = seqs_onehot[np.array(seqs_label) == 0].numpy()
            p_values_per_pwm = []

            for window_size_idx, window_size in enumerate(window_sizes):
                for cluster_idx, pwm in pwms_per_window_size[window_size].items():
                    counts_positive = convolve_pwm_with_sequences(pwm, seqs_positive)
                    counts_negative = convolve_pwm_with_sequences(pwm, seqs_negative)
                    try:
                        chi2, p_value = perform_chi_squared_tests(
                            counts_positive, counts_negative
                        )
                        p_values_per_pwm.append(
                            [
                                pwm,
                                chi2,
                                p_value,
                                cluster_idx,
                                window_size_idx,
                                counts_positive,
                                counts_negative,
                            ]
                        )
                    except ValueError:
                        pass

            corrected_alpha = 0.05 / len(p_values_per_pwm)
            p_values_per_pwm = [_ for _ in p_values_per_pwm if _[2] < corrected_alpha]
            with open(
                f"{self.motif_save_path}/significantly_enriched_motifs.pkl", "wb"
            ) as fp:
                pickle.dump(p_values_per_pwm, fp)

        for pwm, _, p_value, cluster_idx, _, _, _ in p_values_per_pwm:
            self.plot_pwm_logo(pwm, p_value, cluster_idx)
        print(f'Motif plots saved at {self.motif_save_path}.')


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--filename", type=str, choices=["Choi_et_al", "Erlich_et_al", "Gao_et_al", "GCall", "GCfix", "Koch_et_al", "Song_et_al"], required=True)
    parser.add_argument("--threshold", type=str, default='2perc')
    parser.add_argument("--only_visualization", type=bool, default=True)
    args = parser.parse_args()

    motif_analysis = CluMo(args.filename, args.threshold)
    motif_analysis.motif_plot(only_visualization=args.only_visualization)
