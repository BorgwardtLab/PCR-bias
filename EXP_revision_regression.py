from model import CNN_1D_withPE
import pandas as pd
import torch
import os
import warnings
import numpy as np
warnings.filterwarnings("ignore")
import pickle
from torch import nn
import itertools
import random
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    balanced_accuracy_score,
    accuracy_score,
    roc_curve,
    precision_recall_curve
)
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, StratifiedKFold
from utils.training_utils import representation
from torch.utils.data import DataLoader, TensorDataset

def get_ext_data(ext_file_name, threshold):
    seqs = pd.read_pickle(
        "Data/{}/bad_seqs_{}.pkl".format(
            ext_file_name, threshold
        )
    )
    rest_seqs = seqs["rest"]["sequence"]
    bot_seqs = seqs["bottom"]["sequence"]
    rest_qs = seqs["rest"]["eff"]
    bot_qs = seqs["bottom"]["eff"]
    all_seqs = np.hstack([bot_seqs, rest_seqs])
    labels = np.hstack([bot_qs, rest_qs])
    rest_ids = rest_seqs.index.values
    bot_ids = bot_seqs.index.values
    seq_ids = np.hstack([bot_ids, rest_ids])
    X_tar = representation(all_seqs, with_reverse=False)
    X_tar = torch.tensor(X_tar, dtype=torch.float32)
    y_tar = torch.tensor(labels, dtype=torch.float32)

    return X_tar, y_tar, seq_ids
def load_best_model_params(source, threshold):
    best_config = pd.read_pickle(f"CNN/models/{source}_Filtered_{threshold}_best_model_config.pkl")
    return best_config



threshold = '2perc'
filenames = ['GCfix', 'GCall']
criterion = nn.MSELoss()
# # internal validation
for filename in filenames:
    best_params = load_best_model_params(filename, threshold)
    torch.manual_seed(1)
    X_src, y_src, seq_ids = get_ext_data(filename, threshold)
    overall_eff_threshold = np.percentile(y_src.numpy(), 2)
    binary_labels = (y_src.numpy() <= overall_eff_threshold).astype(int)

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=1)

    fold_count = 0
    fold_results = []
    all_preds_df = []
    all_min_pred = float("inf")
    all_max_pred = float("-inf")

    for train_idx, val_idx in cv.split(X_src, binary_labels):
        fold_count += 1

        X_train, X_val = X_src[train_idx], X_src[val_idx]
        y_train, y_val = y_src[train_idx], y_src[val_idx]
        seq_val_ids = seq_ids[val_idx]

        train_dset = TensorDataset(X_train, y_train)
        train_loader = DataLoader(
            train_dset, batch_size=int(best_params["batch_size"]), shuffle=True
        )

        model = CNN_1D_withPE(
            linear_dim=best_params["linear_dim"],
            n_filters=best_params["n_filter"],
            len_filters=best_params["len_filter"],
            normalization=best_params["normalization"],
            norm_first=best_params["norm_first"],
            activation=best_params["activation"],
            AdaPool=best_params["AdaPool"],
            use_PE=True,
            with_reverse=False,
            num_classes=1,
        ).cuda()

        optimizer = torch.optim.Adam(
            model.parameters(),
            lr=best_params["LR"],
            weight_decay=best_params["weight_decay"],
        )

        best_rmse = float("inf")
        best_model = None
        best_val_preds = None
        best_val_true = None

        # Training loop
        for epoch in range(int(150)):
            model.train()
            for data, label in train_loader:
                data, label = data.cuda(), label.cuda()
                optimizer.zero_grad()
                outputs = model(data).squeeze()
                loss = criterion(outputs, label)
                loss.backward()
                optimizer.step()

            # Evaluate on validation set
            model.eval()
            with torch.no_grad():
                val_preds = model(X_val.cuda()).squeeze().cpu().numpy()
                val_true = y_val.numpy()

                rmse = np.sqrt(mean_squared_error(val_true, val_preds))
                if rmse < best_rmse:
                    best_rmse = rmse
                    best_model = model
                    best_val_preds = val_preds
                    best_val_true = val_true

        # After training is done for this fold, compute regression metrics
        fold_r2 = r2_score(best_val_true, best_val_preds)
        fold_mae = mean_absolute_error(best_val_true, best_val_preds)
        fold_rmse = np.sqrt(mean_squared_error(best_val_true, best_val_preds))

        # Track min/max predictions for final probability transform
        all_min_pred = min(all_min_pred, np.min(best_val_preds))
        all_max_pred = max(all_max_pred, np.max(best_val_preds))

        fold_results.append(
            {"fold": fold_count, "r2": fold_r2, "rmse": fold_rmse, "mae": fold_mae}
        )

        df_fold_preds = pd.DataFrame(
            {
                "seq_id": seq_val_ids,
                "true_efficiency": best_val_true,
                "binary_label": (best_val_true <= overall_eff_threshold).astype(
                    int
                ),
                "pred_efficiency": best_val_preds,
                "pred_probability": np.nan,
            }
        )
        df_fold_preds["fold"] = fold_count
        all_preds_df.append(df_fold_preds)

    # Combine fold-level row predictions
    all_preds_df = pd.concat(all_preds_df, axis=0).reset_index(drop=True)
    denom = all_max_pred - all_min_pred if (all_max_pred > all_min_pred) else 1e-9


    def eff_to_prob(eff):
        p = 1.0 - (eff - all_min_pred) / denom
        return np.clip(p, 0.0, 1.0)
    all_preds_df["pred_probability"] = all_preds_df["pred_efficiency"].apply(
        eff_to_prob
    )

    # Save fold regression metrics
    result_df = pd.DataFrame(fold_results)

    # Summaries for regression
    summary_df = pd.DataFrame(
        {
            "r2_mean": [result_df["r2"].mean()],
            "rmse_mean": [result_df["rmse"].mean()],
            "mae_mean": [result_df["mae"].mean()],
        }
    )

    # 5) Visualization: predicted vs true for each fold, and MSE/R^2
    plt.figure()
    plt.scatter(
        all_preds_df["true_efficiency"], all_preds_df["pred_efficiency"], alpha=0.5
    )
    plt.xlabel("True Efficiency")
    plt.ylabel("Predicted Efficiency")
    plt.title(f"{filename}_{threshold} - All Folds Regression Predictions")
    plt.show()

    overall_mse = mean_squared_error(
        all_preds_df["true_efficiency"], all_preds_df["pred_efficiency"]
    )
    overall_r2 = r2_score(
        all_preds_df["true_efficiency"], all_preds_df["pred_efficiency"]
    )
    print(f"Overall MSE: {overall_mse:.4f}, Overall R^2: {overall_r2:.4f}")

    # 6) Compute and visualize ROC / PR for overall prediction using the binary classification logic
    y_all = all_preds_df["binary_label"].values
    y_score = all_preds_df["pred_probability"].values

    # ROC
    fpr, tpr, _ = roc_curve(y_all, y_score)
    roc_auc = roc_auc_score(y_all, y_score)

    # Precision-Recall
    precision, recall, _ = precision_recall_curve(y_all, y_score)
    pr_auc = average_precision_score(y_all, y_score)

    # We'll do per-fold AUROC/AUPRC as well:
    fold_auroc = []
    fold_auprc = []
    for f in range(1, 6):
        mask = all_preds_df["fold"] == f
        y_f = all_preds_df.loc[mask, "binary_label"].values
        s_f = all_preds_df.loc[mask, "pred_probability"].values
        fold_auroc.append(roc_auc_score(y_f, s_f))
        fold_auprc.append(average_precision_score(y_f, s_f))

    # Overall
    plt.figure()
    plt.plot(fpr, tpr, label=f"ROC curve (area = {roc_auc:.3f})")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"{filename}_{threshold} - Overall ROC")
    plt.legend()
    plt.show()

    plt.figure()
    plt.plot(recall, precision, label=f"PR curve (AP = {pr_auc:.3f})")
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title(f"{filename}_{threshold} - Overall Precision-Recall")
    plt.legend()
    plt.show()

    print("Fold-wise AUROC:", fold_auroc)
    print("Fold-wise AUPRC:", fold_auprc)
    print(f"Overall AUROC: {roc_auc:.4f}, Overall AUPRC: {pr_auc:.4f}")
    all_preds_df["seq_id"] = all_preds_df["seq_id"].astype(str).str.lstrip("0")
    all_preds_df.to_csv(
        f"CNN/results_revision/Internal_{filename}_{threshold}_regression_plus_probs.csv",
        index=False,
    )

# external validation
pairs = list(itertools.permutations(filenames, 2))

for source, target in pairs:
    print(f"Running external validation: {source} → {target}")
    best_params = load_best_model_params(source, threshold)

    torch.manual_seed(1)
    X_src, y_src, src_seq_ids = get_ext_data(source, threshold)
    X_tar, y_tar, tar_seq_ids = get_ext_data(target, threshold)

    overall_eff_threshold = np.percentile(y_tar.numpy(), 2)
    binary_labels = (y_tar.numpy() <= overall_eff_threshold).astype(int)

    train_dset = TensorDataset(X_src, y_src)
    train_loader = DataLoader(
        train_dset, batch_size=int(best_params["batch_size"]), shuffle=True
    )

    model = CNN_1D_withPE(
        linear_dim=best_params["linear_dim"],
        n_filters=best_params["n_filter"],
        len_filters=best_params["len_filter"],
        normalization=best_params["normalization"],
        norm_first=best_params["norm_first"],
        activation=best_params["activation"],
        AdaPool=best_params["AdaPool"],
        use_PE=True,
        with_reverse=False,
        num_classes=1,
    ).cuda()

    optimizer = torch.optim.Adam(
        model.parameters(),
        lr=best_params["LR"],
        weight_decay=best_params["weight_decay"],
    )

    best_rmse = float("inf")
    best_model = None
    best_val_preds = None
    best_val_true = None

    # Train on full source
    for epoch in range(150):
        model.train()
        for data, label in train_loader:
            data, label = data.cuda(), label.cuda()
            optimizer.zero_grad()
            outputs = model(data).squeeze()
            loss = torch.nn.functional.mse_loss(outputs, label)
            loss.backward()
            optimizer.step()

    # Evaluate on target
    model.eval()
    with torch.no_grad():
        pred_eff = model(X_tar.cuda()).squeeze().cpu().numpy()
        true_eff = y_tar.numpy()

    # Regression metrics
    rmse = np.sqrt(mean_squared_error(true_eff, pred_eff))
    r2 = r2_score(true_eff, pred_eff)
    mae = mean_absolute_error(true_eff, pred_eff)

    # Store predictions
    df_preds = pd.DataFrame({
        "seq_id": tar_seq_ids,
        "true_efficiency": true_eff,
        "binary_label": binary_labels,
        "pred_efficiency": pred_eff,
        "pred_probability": np.nan,
    })

    # Map predictions to probabilities (same transformation logic)
    min_pred = np.min(pred_eff)
    max_pred = np.max(pred_eff)
    denom = max_pred - min_pred if max_pred > min_pred else 1e-9
    def eff_to_prob(eff):
        return np.clip(1.0 - (eff - min_pred) / denom, 0.0, 1.0)

    df_preds["pred_probability"] = df_preds["pred_efficiency"].apply(eff_to_prob)

    # Plots
    plt.figure()
    plt.scatter(true_eff, pred_eff, alpha=0.5)
    plt.xlabel("True Efficiency")
    plt.ylabel("Predicted Efficiency")
    plt.title(f"{source}→{target}_{threshold} - External Prediction")
    plt.show()

    # ROC / PR
    y_true = df_preds["binary_label"].values
    y_score = df_preds["pred_probability"].values

    fpr, tpr, _ = roc_curve(y_true, y_score)
    precision, recall, _ = precision_recall_curve(y_true, y_score)
    roc_auc = roc_auc_score(y_true, y_score)
    pr_auc = average_precision_score(y_true, y_score)

    plt.figure()
    plt.plot(fpr, tpr, label=f"ROC (AUC = {roc_auc:.3f})")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"{source}→{target} - External ROC")
    plt.legend()
    plt.show()

    plt.figure()
    plt.plot(recall, precision, label=f"PR (AP = {pr_auc:.3f})")
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title(f"{source}→{target} - External PR")
    plt.legend()
    plt.show()

    # Print metrics
    print(f"[{source}→{target}] R2: {r2:.4f}, RMSE: {rmse:.4f}, MAE: {mae:.4f}")
    print(f"AUROC: {roc_auc:.4f}, AUPRC: {pr_auc:.4f}")

    # Save
    df_preds["seq_id"] = df_preds["seq_id"].astype(str).str.lstrip("0")
    df_preds.to_csv(
        f"CNN/results_revision/External_{source}2{target}_{threshold}_regression_plus_probs.csv",
        index=False,
    )