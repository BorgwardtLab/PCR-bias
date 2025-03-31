import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader, TensorDataset
import pytorch_lightning as pl

all_datasets = [
    "Choi_et_al",
    "Erlich_et_al",
    "Gao_et_al",
    "GCall",
    "GCfix",
    "Koch_et_al",
    "Song_et_al",
]

def get_data(ext_file_name: str, threshold: str):
    seqs = pd.read_pickle("Data/{}/bad_seqs_{}.pkl".format(ext_file_name, threshold))
    rest_seqs = seqs["rest"]["sequence"]
    bot_seqs = seqs["bottom"]["sequence"]
    seqs = np.hstack([bot_seqs, rest_seqs])
    labels = [1] * len(bot_seqs) + [0] * len(rest_seqs)
    X_tar = representation(seqs, with_reverse=False)
    X_tar = torch.tensor(X_tar)
    y_tar = torch.tensor(labels)
    return X_tar, y_tar


def reverse_complement(dna):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join([complement[base] for base in dna[::-1] if base in "ACGT"])


def get_DNA_seq_array_reverse_complement(seq):
    alpha = "ACGT"
    row = len(seq)
    new_array = np.zeros((row, 4))
    reverse_array = np.zeros((row, 4))
    rev_seq = reverse_complement(seq)

    for i, (val, rev_val) in enumerate(zip(seq, rev_seq)):
        if val in alpha:
            index = alpha.index(val)
            new_array[i][index] = 1
        if rev_val in alpha:
            rev_index = alpha.index(rev_val)
            reverse_array[i][rev_index] = 1
    new_array = np.hstack((new_array, reverse_array))
    return new_array


def get_DNA_seq_array(seq):
    alpha = "ACGT"
    row = len(seq)
    new_array = np.zeros((row, 4))

    for i, val in enumerate(seq):
        index = alpha.index(val)
        new_array[i][index] = 1

    return new_array


def representation(X, with_reverse):
    reformed_seqs = []
    for i in range(len(X)):
        if with_reverse:
            seq = get_DNA_seq_array_reverse_complement(X[i])
        else:
            seq = get_DNA_seq_array(X[i])
        reformed_seqs.append(seq)
    reformed_seqs = np.array(reformed_seqs)
    return reformed_seqs

class DNADataModule(pl.LightningDataModule):
    def __init__(self, X_train, y_train, X_val, y_val, batch_size=64):
        super().__init__()
        self.X_train = X_train
        self.y_train = y_train
        self.X_val = X_val
        self.y_val = y_val
        self.batch_size = batch_size

    def setup(self, stage=None):
        # Just store PyTorch TensorDatasets
        self.train_ds = TensorDataset(self.X_train, self.y_train)
        self.val_ds   = TensorDataset(self.X_val,   self.y_val)

    def train_dataloader(self):
        return DataLoader(self.train_ds, batch_size=self.batch_size, shuffle=True)

    def val_dataloader(self):
        return DataLoader(self.val_ds, batch_size=self.batch_size, shuffle=False)