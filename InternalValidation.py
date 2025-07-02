import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import torch
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import average_precision_score, roc_auc_score
import pytorch_lightning as pl
from pytorch_lightning.callbacks import ModelCheckpoint
from model import CNN
from utils.training_utils import get_data, DNADataModule


def parse_arguments():
    parser = argparse.ArgumentParser(description="1D CNN Internal Validation")
    parser.add_argument("--filename", type=str, choices=["Choi_et_al", "Erlich_et_al", "Gao_et_al", "GCall", "GCfix", "Koch_et_al", "Song_et_al"], required=True)
    parser.add_argument("--threshold", type=str, default='2perc')
    return parser.parse_args()


def load_data(filename, threshold):
    torch.manual_seed(1)
    return get_data(filename, threshold)


def load_best_model_params(source, threshold):
    best_config = pd.read_pickle(f"CNN/models/{source}_Filtered_{threshold}_best_model_config.pkl")
    return best_config


def train_and_evaluate_model(
    X_src, y_src, params, fold_count, filename, threshold, train_idx, val_idx
):
    X_train = X_src[train_idx]
    y_train = y_src[train_idx]
    X_val = X_src[val_idx]
    y_val = y_src[val_idx]

    data_module = DNADataModule(
        X_train=X_train,
        y_train=y_train,
        X_val=X_val,
        y_val=y_val,
        batch_size=params["batch_size"],
    )
    model = CNN(hparams=params)
    checkpoint_callback = ModelCheckpoint(
        monitor="val_auprc",
        mode="max",
        save_top_k=1,
        dirpath=f"CNN/checkpoints/internal/{filename}/{threshold}/",
        filename=f"best_model_fold_{fold_count}",
    )

    trainer = pl.Trainer(
        max_epochs=params["n_epoch"],
        gpus=-1 if torch.cuda.is_available() else 0,
        callbacks=[checkpoint_callback],
        enable_progress_bar=False,
    )

    trainer.fit(model, data_module)

    best_model_path = checkpoint_callback.best_model_path
    best_model = CNN.load_from_checkpoint(best_model_path)

    best_model.eval()
    with torch.no_grad():
        preds = torch.sigmoid(best_model(X_val.float()))[:, 1].cpu().numpy()
        labels = y_val.numpy()

    auprc = average_precision_score(labels, preds)
    auroc = roc_auc_score(labels, preds)

    return {
        "fold": fold_count,
        "test_auprc": auprc,
        "test_auroc": auroc,
        "y_test": list(labels),
        "pred_test": list(preds),
    }


def save_metrics(metrics_summary, filename, threshold):
    result_path = Path(f"CNN/results/internal/{filename}/{threshold}/")
    result_path.mkdir(parents=True, exist_ok=True)

    df_metrics = pd.DataFrame(metrics_summary).drop(columns=["y_test", "pred_test"])
    df_metrics.to_csv(f"{result_path}/1DCNN_PE_per_fold.csv", index=False)

    all_preds = np.concatenate([m["pred_test"] for m in metrics_summary])
    all_labels = np.concatenate([m["y_test"] for m in metrics_summary])

    overall_auroc = roc_auc_score(all_labels, all_preds)
    overall_auprc = average_precision_score(all_labels, all_preds)

    test_aurocs = [m["test_auroc"] for m in metrics_summary]
    test_auprcs = [m["test_auprc"] for m in metrics_summary]

    summary_stats = {
        "mean_test_auroc": np.mean(test_aurocs),
        "std_test_auroc": np.std(test_aurocs),
        "mean_test_auprc": np.mean(test_auprcs),
        "std_test_auprc": np.std(test_auprcs),
        "pooled_test_auroc": overall_auroc,
        "pooled_test_auprc": overall_auprc,
    }

    df_summary = pd.DataFrame([summary_stats])
    df_summary.to_csv(f"{result_path}/1DCNN_PE_summary.csv", index=False)


def main():
    args = parse_arguments()
    filename, threshold = args.filename, args.threshold
    X_src, y_src = load_data(filename, threshold)
    best_config = load_best_model_params(filename, threshold)
    cv = StratifiedKFold(n_splits=5, shuffle=False)
    metrics_summary = []
    for fold_count, (train_idx, val_idx) in enumerate(cv.split(X_src, y_src), start=1):
        fold_metrics = train_and_evaluate_model(
            X_src, y_src, best_config, fold_count, filename, threshold, train_idx, val_idx
        )
        metrics_summary.append(fold_metrics)
    save_metrics(metrics_summary, filename, threshold)
    print(f'Results saved at CNN/results/internal/{filename}/{threshold}/')

if __name__ == "__main__":
    main()
