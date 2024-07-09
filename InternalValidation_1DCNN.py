import argparse
from pathlib import Path
import os
import pandas as pd
import numpy as np
import torch
from sklearn.model_selection import StratifiedKFold, train_test_split, ParameterSampler
from sklearn.metrics import average_precision_score, roc_auc_score
import pytorch_lightning as pl
from pytorch_lightning.callbacks import ModelCheckpoint
from model import CNN
from utils.training_utils import CNN_hyper_params, get_data, DNADataModule


def parse_arguments():
    parser = argparse.ArgumentParser(description="1D CNN Internal Validation")
    parser.add_argument("--filename", type=str, required=True)
    parser.add_argument("--threshold", type=str, required=True)
    return parser.parse_args()


def load_data(filename, threshold):
    torch.manual_seed(1)
    return get_data(filename, threshold)

def train_and_evaluate_model(X_src, y_src, params, fold_count, model_id, filename, threshold, train_idx, test_idx):
    X_train, X_test = X_src[train_idx], X_src[test_idx]
    y_train, y_test = y_src[train_idx], y_src[test_idx]

    X_train_inner, X_val, y_train_inner, y_val = train_test_split(
        X_train, y_train, test_size=0.1 / 0.8, stratify=y_train
    )

    model = CNN(hparams=params)
    data_module = DNADataModule(
        X_train=X_train_inner,
        y_train=y_train_inner,
        X_val=X_val,
        y_val=y_val,
        X_test=X_test,
        y_test=y_test,
        batch_size=params["batch_size"],
    )

    checkpoint_callback = ModelCheckpoint(
        monitor="val_auprc",
        mode="max",
        save_top_k=1,
        dirpath=f"CNN/model/internal/{filename}/{threshold}/",
        filename=f"best_model_fold_{fold_count}_model_id_{model_id}",
    )

    trainer = pl.Trainer(
        max_epochs=100,
        gpus=-1 if torch.cuda.is_available() else 0,
        callbacks=[checkpoint_callback],
        progress_bar_refresh_rate=20,
    )

    trainer.fit(model, data_module)

    best_model_path = checkpoint_callback.best_model_path
    best_model = CNN.load_from_checkpoint(best_model_path)

    pred_test = torch.sigmoid(best_model(X_test.float()))[:, 1]
    test_auprc = average_precision_score(y_test.detach().numpy(), pred_test.detach().numpy())
    test_auroc = roc_auc_score(y_test.detach().numpy(), pred_test.detach().numpy())

    pred_val = torch.sigmoid(best_model(X_val.float()))[:, 1]
    val_auprc = average_precision_score(y_val.detach().numpy(), pred_val.detach().numpy())
    val_auroc = roc_auc_score(y_val.detach().numpy(), pred_val.detach().numpy())

    return {
        "model id": model_id,
        "fold": fold_count,
        **params,
        "val_auprc": val_auprc,
        "val_auroc": val_auroc,
        "test_auprc": test_auroc,
        "test_auroc": test_auprc,
    }


def remove_non_best_models(metrics_summary, best_model_id, filename, threshold):
    dirpath = f"CNN/model/internal/{filename}/{threshold}/"
    for model in metrics_summary:
        if model["model id"] != best_model_id:
            for fold_count in range(5):
                model_path = Path(dirpath) / f"best_model_fold_{fold_count}_model_id_{model['model id']}.ckpt"
                if model_path.exists():
                    model_path.unlink()
    for fold_count in range(5):
        model_old = Path(dirpath) / f"best_model_fold_{fold_count}_model_id_{best_model_id}.ckpt"
        if model_old.exists():
            model_new = Path(dirpath) / f"best_model_fold_{fold_count}.ckpt"
            os.rename(model_old, model_new)


def save_metrics(metrics_summary, filename, threshold):
    metrics_df = pd.DataFrame(metrics_summary)
    result_path = Path(f'CNN/results/internal/{filename}/{threshold}/')
    result_path.mkdir(parents=True, exist_ok=True)
    metrics_df.to_csv(result_path / '1DCNN_PE.csv', index=False)


def main():
    args = parse_arguments()
    filename, threshold = args.filename, args.threshold

    X_src, y_src = load_data(filename, threshold)
    param_sampler = ParameterSampler(CNN_hyper_params, n_iter=50)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=1)

    metrics_summary = []
    avg_val_auprc = {}

    for model_id, params in enumerate(param_sampler, start=1):
        if params.get('use_class_weight') and threshold == '2perc':
            params['threshold'] = 2

        val_auprcs = []

        for fold_count, (train_idx, test_idx) in enumerate(cv.split(X_src, y_src)):
            metrics = train_and_evaluate_model(X_src, y_src, params, fold_count, model_id, filename, threshold, train_idx, test_idx)
            metrics_summary.append(metrics)
            val_auprcs.append(metrics["val_auprc"])

        avg_val_auprc[model_id] = np.mean(val_auprcs)

    best_model_id = max(avg_val_auprc, key=avg_val_auprc.get)
    remove_non_best_models(metrics_summary, best_model_id, filename, threshold)
    save_metrics(metrics_summary, filename, threshold)


if __name__ == "__main__":
    main()
