"""Thin wrappers around the CPA package."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd

from .cpa_prep import CONTROL_GROUP, CPA_DOSE_KEY


def ensure_cpa_compat() -> None:
    compat = Path(__file__).resolve().parents[1] / "cpa_compat"
    if compat.exists() and str(compat) not in sys.path:
        sys.path.insert(0, str(compat))


def setup_cpa_anndata(adata, use_rdkit_embeddings: bool = True):
    ensure_cpa_compat()
    import cpa

    dosage_key = CPA_DOSE_KEY if CPA_DOSE_KEY in adata.obs else "log_dose"
    setup_kwargs = {
        "perturbation_key": "condition_ID",
        "control_group": CONTROL_GROUP,
        "dosage_key": dosage_key,
        "is_count_data": False,
        "categorical_covariate_keys": ["cell_type"],
    }
    if use_rdkit_embeddings:
        setup_kwargs["smiles_key"] = "smiles_rdkit"
    else:
        cpa.CPA.pert_smiles_map = None
    cpa.CPA.setup_anndata(adata, **setup_kwargs)
    return adata


def default_hyperparams() -> dict:
    return {
        "recon_loss": "gauss",
        "doser_type": "logsigm",
        "n_latent": 128,
        "n_hidden_encoder": 512,
        "n_layers_encoder": 3,
        "n_hidden_decoder": 512,
        "n_layers_decoder": 3,
        "n_hidden_doser": 128,
        "n_layers_doser": 2,
        "dropout_rate_encoder": 0.25,
        "dropout_rate_decoder": 0.25,
        "variational": False,
    }


def train_cpa(
    h5ad: str | Path,
    out_dir: str | Path,
    max_epochs: int,
    batch_size: int,
    use_gpu: bool | str | int = True,
    check_val_every_n_epoch: int = 5,
    early_stopping_patience: int = 10,
    hyperparams: dict | None = None,
    plan_kwargs: dict | None = None,
    use_rdkit_embeddings: bool = True,
):
    ensure_cpa_compat()
    import scanpy as sc
    import cpa

    adata = sc.read_h5ad(h5ad)
    setup_cpa_anndata(adata, use_rdkit_embeddings=use_rdkit_embeddings)
    params = default_hyperparams()
    if hyperparams:
        params.update(hyperparams)
    model = cpa.CPA(
        adata,
        split_key="split",
        train_split="train",
        valid_split="test",
        test_split="ood",
        use_rdkit_embeddings=use_rdkit_embeddings,
        **params,
    )
    model.train(
        max_epochs=max_epochs,
        use_gpu=use_gpu,
        batch_size=batch_size,
        check_val_every_n_epoch=check_val_every_n_epoch,
        early_stopping_patience=early_stopping_patience,
        plan_kwargs=plan_kwargs,
        save_path=str(out_dir),
    )
    return model


def load_cpa(
    model_dir: str | Path,
    h5ad: str | Path,
    use_gpu: bool | str | int = True,
    use_rdkit_embeddings: bool = True,
):
    ensure_cpa_compat()
    import scanpy as sc
    import cpa

    adata = sc.read_h5ad(h5ad)
    setup_cpa_anndata(adata, use_rdkit_embeddings=use_rdkit_embeddings)
    return cpa.CPA.load(str(model_dir), adata=adata, use_gpu=use_gpu), adata


def predict_query_rows(
    model_dir: str | Path,
    h5ad: str | Path,
    out_h5ad: str | Path,
    out_target_long: str | Path | None = None,
    targets: str | Path | None = None,
    response_source: str = "final_cpa_unknown_prediction",
    batch_size: int = 512,
    use_gpu: bool | str | int = True,
    use_rdkit_embeddings: bool = True,
) -> None:
    import anndata as ad

    model, adata = load_cpa(model_dir, h5ad, use_gpu=use_gpu, use_rdkit_embeddings=use_rdkit_embeddings)
    query_idx = np.where(adata.obs["split"].astype(str).eq("query"))[0]
    if len(query_idx) == 0:
        raise ValueError("No rows with split == 'query' found in CPA AnnData")
    pred = model.custom_predict(adata=adata, indices=query_idx, batch_size=batch_size, n_samples=1)["latent_x_pred"]
    pred.var = adata.var.copy()
    pred.obs = adata.obs.iloc[query_idx].copy()
    out_h5ad = Path(out_h5ad)
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    pred.write_h5ad(out_h5ad, compression="gzip")
    if out_target_long is not None:
        target_long = predictions_to_target_long(pred, targets=targets, response_source=response_source)
        Path(out_target_long).parent.mkdir(parents=True, exist_ok=True)
        target_long.to_csv(out_target_long, sep="\t", index=False)


def predictions_to_target_long(
    pred,
    targets: str | Path | None = None,
    response_source: str = "final_cpa_unknown_prediction",
) -> pd.DataFrame:
    if targets is not None:
        target_df = pd.read_csv(targets, sep="\t")
        genes = target_df["gene_name"].dropna().astype(str).unique().tolist()
    else:
        genes = pred.var_names.astype(str).tolist()
    genes = [g for g in genes if g in set(pred.var_names.astype(str))]
    if not genes:
        raise ValueError("No requested target genes are present in the CPA prediction matrix")
    sub = pred[:, genes]
    x = np.asarray(sub.X, dtype=float)
    frames = []
    meta_cols = ["pert_id", "pert_iname", "condition_ID", "cell_id", "pert_dose", "dose_um", "log_dose", CPA_DOSE_KEY]
    meta_cols = [c for c in meta_cols if c in sub.obs.columns]
    for j, gene in enumerate(sub.var_names.astype(str)):
        frame = sub.obs[meta_cols].copy()
        frame["gene_name"] = gene
        frame["mean_z"] = x[:, j]
        frame["median_z"] = x[:, j]
        frame["n_signatures"] = 1
        frame["response_source"] = response_source
        frames.append(frame)
    out = pd.concat(frames, ignore_index=True)
    return out
