import argparse
from pathlib import Path
import pandas as pd
import torch
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import train_test_split
import pytorch_lightning as pl
from pytorch_lightning.callbacks import ModelCheckpoint
from model import CNN
from utils.training_utils import get_data, DNADataModule, all_datasets
def parse_arguments():
    parser = argparse.ArgumentParser(description="1D CNN External Validation")
    parser.add_argument("--filename", type=str, required=True)
    parser.add_argument("--threshold", type=str, required=True)
    return parser.parse_args()

def load_best_model_params(source, threshold):
    metrics_path = Path(f'CNN/results/internal/{source}/{threshold}/1DCNN_PE.csv')
    metrics = pd.read_csv(metrics_path)
    best_model_id = metrics.groupby('model id')['val_auprc'].mean().idxmax()
    return metrics.loc[metrics['model id'] == best_model_id].iloc[0].to_dict()

def train_model(source, threshold, best_model_params, X_train, y_train, X_val, y_val):
    model = CNN(hparams=best_model_params)
    data_module = DNADataModule(
        X_train=X_train,
        y_train=y_train,
        X_val=X_val,
        y_val=y_val,
        batch_size=best_model_params["batch_size"],
    )

    checkpoint_callback = ModelCheckpoint(
        monitor="val_auprc",
        mode="max",
        save_top_k=1,
        dirpath=f"CNN/model/external/{source}/{threshold}/",
        filename=f"best_model",
    )

    trainer = pl.Trainer(
        max_epochs=100,
        gpus=-1 if torch.cuda.is_available() else 0,
        callbacks=[checkpoint_callback],
        progress_bar_refresh_rate=20,
    )

    trainer.fit(model, data_module)
    return checkpoint_callback.best_model_path

def evaluate_model_on_targets(model, source, threshold):
    metrics_summary = []
    for target in all_datasets:
        if source != target:
            X_tar, y_tar = get_data(target, threshold)
            with torch.no_grad():
                pred_tar = torch.sigmoid(model(X_tar.float()))[:, 1]
            tar_auprc = average_precision_score(y_tar.numpy(), pred_tar.numpy())
            tar_auroc = roc_auc_score(y_tar.numpy(), pred_tar.numpy())
            metrics_summary.append(
                {
                    "target": target,
                    "test_auprc": tar_auprc,
                    "test_auroc": tar_auroc,
                }
            )
    return pd.DataFrame(metrics_summary)

def save_metrics(metrics_df, source, threshold):
    result_path = Path(f'CNN/results/external/{source}/{threshold}/')
    result_path.mkdir(parents=True, exist_ok=True)
    metrics_df.to_csv(result_path / '1DCNN_PE.csv')

def main():
    args = parse_arguments()
    torch.manual_seed(1)
    X_src, y_src = get_data(args.filename, args.threshold)
    best_model_params = load_best_model_params(args.filename, args.threshold)
    X_train, X_val, y_train, y_val = train_test_split(X_src, y_src, test_size=0.2, stratify=y_src)

    best_model_path = train_model(args.filename, args.threshold, best_model_params, X_train, y_train, X_val, y_val)
    best_model = CNN.load_from_checkpoint(best_model_path)
    metrics_df = evaluate_model_on_targets(best_model, args.filename, args.threshold)
    save_metrics(metrics_df, args.filename, args.threshold)

if __name__ == "__main__":
    main()
