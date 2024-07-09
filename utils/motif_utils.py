import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logomaker
from scipy.stats import entropy, chi2_contingency
import re
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score
from sklearn.cluster import KMeans

def one_hot_encode(seq):
    allowed = set("ACTG")
    if not set(seq).issubset(allowed):
        invalid = set(seq) - allowed
        raise ValueError(
            f"Sequence contains chars not in allowed DNA alphabet (ACGT): {invalid}"
        )

    nuc_d = {
        "A": [1.0, 0.0, 0.0, 0.0],
        "C": [0.0, 1.0, 0.0, 0.0],
        "G": [0.0, 0.0, 1.0, 0.0],
        "T": [0.0, 0.0, 0.0, 1.0],
    }
    vec = np.array([nuc_d[x] for x in seq])
    return vec

def one_hot_to_DNA(one_hot_encoded_seq):
    alpha = "ACGT"
    dna_seq = ""
    for row in one_hot_encoded_seq:
        if 1 in row:
            index = np.argmax(row)
            dna_seq += alpha[index]
    return dna_seq
def tsne_kmeans_search_silhouette(distance_mat, sample_weight=None):
    k_values = range(2, 10)
    perplexities = [2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20]
    best_k = None
    best_perplexity = None
    best_embeddings = None
    best_score = -1
    for perplexity in perplexities:
        tsne = TSNE(
            n_components=2,
            perplexity=perplexity,
            metric="precomputed",
            random_state=42,
            init="random",
        )
        embeddings = tsne.fit_transform(distance_mat)
        for k in k_values:
            kmeans = KMeans(n_clusters=k, random_state=42)
            labels = kmeans.fit_predict(embeddings, sample_weight=sample_weight)
            if len(set(labels)) > 1:
                score = silhouette_score(embeddings, labels)
                if score > best_score:
                    best_k = k
                    best_score = score
                    best_embeddings = embeddings
                    best_perplexity = perplexity

    best_kmeans = KMeans(n_clusters=best_k, random_state=42)
    return best_kmeans.fit_predict(best_embeddings), best_embeddings


def find_most_significant_kmer_in_sequence(
    attribution_matrix,
    window_size=4,
    stride=1,
):
    best_score = -np.inf
    best_motif_start = None
    attribution_threshold = np.max(attribution_matrix) / 2
    average_score = []
    for start_pos in range(0, attribution_matrix.shape[0] - window_size + 1, stride):
       
        window = attribution_matrix[start_pos : start_pos + window_size, :]
        if len(window) == window_size:
            sum_score = np.sum(window)
            average_score.append(sum_score)
            if sum_score > best_score:
                best_score = sum_score
                best_motif_start = start_pos
    if best_score > 2 * np.mean(average_score):
        return best_motif_start, best_score
    else:
        return None, None

def one_hot_to_DNA(one_hot_encoded_seq):
    alpha = "ACGT"
    dna_seq = ""
    for row in one_hot_encoded_seq:
        if 1 in row:
            index = np.argmax(row)
            dna_seq += alpha[index]
    return dna_seq


def hamming_distance(x, y):
    return sum(c1 != c2 for c1, c2 in zip(x, y))

def pairwise_distance(strings):
    n = len(strings)
    distance_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            distance = hamming_distance(strings[i], strings[j])
            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance
    return distance_matrix

def count_motif_in_segments(motif, sequences, region):
    motif_pattern = re.compile(motif)
    counts = [0, 0]
    for seq in sequences:
        segment = seq[region[0] : region[1]]
        if motif_pattern.search(segment):
            counts[0] += 1
        else:
            counts[1] += 1
    return counts
def create_pwm_for_cluster(cluster_df):
    nucleotide_to_num = {"A": 0, "C": 1, "G": 2, "T": 3}
    pwm = np.zeros((4, len(cluster_df["sequence"].iloc[0])))

    # Calculate weighted counts
    for _, row in cluster_df.iterrows():
        for i, nucleotide in enumerate(row["sequence"]):
            if i < len(row["sequence"]):
                pwm[nucleotide_to_num[nucleotide], i] += row["occurrence"]
            else:
                pass
    pwm /= pwm.sum(axis=0)
    return pwm


def perform_chi_squared_tests(counts_positive, counts_negative):
    contingency_table = np.array([counts_positive, counts_negative])
    chi2, p = chi2_contingency(contingency_table)[:2]
    return chi2, p

def convolve_pwm_with_sequences(pwm, sequences):
    pwm_length = pwm.shape[1]
    convolved_scores = 0

    for seq in sequences:
        sequence_length = seq.shape[0]
        scores = []

        # Convolution step
        for i in range(sequence_length - pwm_length + 1):
            window = seq[i : i + pwm_length, :]
            score = np.sum(pwm * window.T) / pwm_length
            scores.append(score)

        # Max pooling step
        max_score = np.max(scores) if scores else 0
        if max_score > 0.5:
            convolved_scores += 1
    counts = np.array(
        [np.ceil(convolved_scores), len(sequences) - np.ceil(convolved_scores)]
    )
    return counts
