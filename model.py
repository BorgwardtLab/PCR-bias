import numpy as np
import torch
import torchmetrics
from torch import nn
import pytorch_lightning as pl

activation_dict = {
    "relu": nn.ReLU(),
    "elu": nn.ELU(),
    "leakyrelu": nn.LeakyReLU(),
    "tanh": nn.Tanh(),
    "sigmoid": nn.Sigmoid(),
}

normalization_dict = {"batchnorm": nn.BatchNorm1d, "instancenorm": nn.InstanceNorm1d}

pooling_dic = {"Avg": nn.AdaptiveAvgPool1d(1), "Max": nn.AdaptiveMaxPool1d(1)}


class CNN_1D_withPE(nn.Module):
    def __init__(
        self,
        linear_dim,
        n_filters,
        len_filters,
        normalization,
        norm_first,
        activation,
        AdaPool,
        use_PE=False,
        with_reverse=False,
        num_classes=2,
    ):
        super().__init__()

        if with_reverse:
            input_dim = 8
        else:
            input_dim = 4

        activation = activation_dict[activation]

        # Encoder Layers
        self.encoder = nn.ModuleList(
            [
                nn.Linear(input_dim, linear_dim),
                nn.Conv1d(linear_dim, n_filters, len_filters),
                activation,
            ]
        )

        if normalization is not None:
            normalization = normalization_dict[normalization]
            if norm_first:
                self.encoder.insert(
                    2, normalization(n_filters)
                )  # After first activation
            else:
                self.encoder.append(normalization(n_filters))  # After convolution

        # Pooling Layer
        if AdaPool == "Avg":
            self.pool = nn.AdaptiveAvgPool1d(1)
        elif AdaPool == "Max":
            self.pool = nn.AdaptiveMaxPool1d(1)

        # Classifier Layer
        self.classifier = nn.Linear(n_filters, num_classes)

        # Positional Encoding Dimension
        self.use_PE = use_PE
        self.linear_dim = linear_dim

    def generate_positional_encodings(self, sequence_length, x):
        position = np.arange(sequence_length)[:, np.newaxis]
        div_term = np.exp(
            np.arange(0, self.linear_dim, 2) * -(np.log(10000.0) / self.linear_dim)
        )
        positional_encodings = np.zeros((sequence_length, self.linear_dim))
        positional_encodings[:, 0::2] = np.sin(position * div_term)
        positional_encodings[:, 1::2] = np.cos(position * div_term)

        positional_encodings = torch.tensor(positional_encodings)
        positional_encodings = (
            torch.swapaxes(positional_encodings, 0, 1).unsqueeze(0).float()
        )
        positional_encodings = positional_encodings.repeat(x.size(0), 1, 1)
        return positional_encodings

    def forward(self, x):
        # Linear projection
        x = self.encoder[0](x)
        x = torch.swapaxes(x, 1, 2)
        # Add positional encodings if use_PE is True
        if self.use_PE:
            pe = self.generate_positional_encodings(x.shape[2], x).float().to(x.device)
            # pe = self.generate_positional_encodings(x.shape[2], x).float()
            x = x + pe

        # Encoding through the rest of the layers
        for layer in self.encoder[1:]:
            x = layer(x)

        # Pooling
        if self.pool is not None:
            x = self.pool(x)

        # Squeeze and Classifying
        x = torch.squeeze(x)
        x = self.classifier(x)

        return x


class CNN(pl.LightningModule):
    def __init__(self, hparams):
        super().__init__()
        self.save_hyperparameters(hparams)
        self.model = CNN_1D_withPE(
            linear_dim=self.hparams["linear_dim"],
            n_filters=self.hparams["n_filter"],
            len_filters=self.hparams["len_filter"],
            normalization=self.hparams["normalization"],
            norm_first=self.hparams["norm_first"],
            activation= self.hparams["activation"],
            AdaPool=self.hparams["AdaPool"],
            use_PE=self.hparams["use_PE"],
            num_classes=2,
        )

        self.criterion = nn.CrossEntropyLoss()

        # Validation metrics
        self.val_auprc = torchmetrics.AveragePrecision(num_classes=2)
        self.val_auroc = torchmetrics.AUROC(num_classes=2)

    def forward(self, x):
        return self.model(x)

    def training_step(self, batch, batch_idx):
        x, y = batch
        logits = self(x.float())
        loss = self.criterion(logits, y)
        self.log("train_loss", loss, on_step=False, on_epoch=True, prog_bar=True)
        return loss

    def validation_step(self, batch, batch_idx):
        x, y = batch
        logits = self(x.float())
        loss = self.criterion(logits, y)
        self.val_auprc(logits.softmax(dim=-1), y)
        self.val_auroc(logits.softmax(dim=-1), y)
        return {"val_loss": loss, "probs": logits.softmax(dim=-1), "targets": y}

    def validation_epoch_end(self, outputs):
        auprc = self.val_auprc.compute()
        auroc = self.val_auroc.compute()
        self.log("val_auprc", auprc, on_epoch=True, prog_bar=True)
        self.log("val_auroc", auroc, on_epoch=True, prog_bar=True)
        self.val_auprc.reset()
        self.val_auroc.reset()

    def configure_optimizers(self):
        return torch.optim.Adam(
            self.parameters(), lr=self.hparams["LR"], weight_decay=self.hparams["weight_decay"]
        )
