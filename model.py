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

normalization_dict = {
    "batchnorm": nn.BatchNorm1d,
    "instancenorm": nn.InstanceNorm1d
}

pooling_dic = {
    "Avg": nn.AdaptiveAvgPool1d(1),
    "Max": nn.AdaptiveMaxPool1d(1)
}

class CNN_1D_withPE(nn.Module):

    def __init__(
        self,
        number_layers: int,
        linear_dim: int,
        n_filters: int,
        len_filters: int,
        normalization: str,
        activation: str,
        AdaPool: str,
        use_PE: bool,
        num_classes: int,
    ):
        super().__init__()
        input_dim = 4
        activation = activation_dict[activation]
        normalization = normalization_dict[normalization]

        self.linear_dim = linear_dim
        self.use_PE = use_PE
        self.linear_layer = nn.Linear(input_dim, self.linear_dim)
        self.encoder = nn.ModuleList()
        for i in range(number_layers):
            if i == 0:
                in_dim = linear_dim
            else:
                in_dim = n_filters
            self.encoder.append(nn.Conv1d(in_dim, n_filters, len_filters))
            self.encoder.append(activation)
            self.encoder.append(normalization(n_filters))

        self.pool = pooling_dic[AdaPool]
        self.classifier = nn.Linear(n_filters, num_classes)


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
        # Linear embedding
        x = torch.swapaxes(self.linear_layer(x), 1, 2)
        # Add positional encodings if use_PE set to True
        if self.use_PE:
            pe = self.generate_positional_encodings(x.shape[2], x).float().to(x.device)
            x = x + pe
        # Encoding through the rest of the layers
        for layer in self.encoder:
            x = layer(x)
        # Pooling
        if self.pool is not None:
            x = self.pool(x)
        # Classifier
        x = self.classifier(torch.squeeze(x))
        return x


class CNN(pl.LightningModule):
    def __init__(self, hparams):
        super().__init__()
        self.save_hyperparameters(hparams)
        self.model = CNN_1D_withPE(
            number_layers=self.hparams['num_layers'],
            linear_dim=self.hparams['linear_dim'],
            n_filters=self.hparams['n_filters'],
            len_filters=self.hparams['len_filters'],
            normalization="batchnorm",
            activation="relu",
            AdaPool=self.hparams['AdaPool'],
            use_PE=self.hparams['use_PE'],
            num_classes=2
        )
        if self.hparams['use_class_weight']:
            self.criterion = nn.CrossEntropyLoss(
                weight=torch.tensor([self.hparams['threshold'], 100 - self.hparams['threshold']]).float().to(
                    self.device))
        else:
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
        self.log('train_loss', loss, on_step=False, on_epoch=True, prog_bar=True)
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
        self.log('val_auprc', auprc, on_epoch=True, prog_bar=True)
        self.log('val_auroc', auroc, on_epoch=True, prog_bar=True)
        self.val_auprc.reset()
        self.val_auroc.reset()

    def configure_optimizers(self):
        return torch.optim.Adam(self.parameters(), lr=self.hparams['LR'], weight_decay=self.hparams['wd'])


