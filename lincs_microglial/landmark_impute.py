"""Linear landmark-to-all-gene imputation for CPA fallback outputs."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from .lincs_io import DEFAULT_GCTX, DEFAULT_GENE_INFO, read_gene_info, read_gctx_block, read_sig_info


def fit_ridge_landmark_to_all(
    gctx: str | Path = DEFAULT_GCTX,
    sig_info: str | Path = "data/raw/lincs_gse92742/GSE92742_Broad_LINCS_sig_info.txt.gz",
    gene_info: str | Path = DEFAULT_GENE_INFO,
    out_npz: str | Path = "data/processed/cpa/landmark_to_all_ridge.npz",
    max_signatures: int = 50000,
    alpha: float = 100.0,
    seed: int = 1,
    chunk_size: int = 1024,
) -> None:
    sig = read_sig_info(sig_info)
    sig = sig[sig["pert_type"].isin(["trt_cp", "ctl_vehicle", "ctl_vehicle.cns"])].copy()
    if len(sig) > max_signatures:
        sig = sig.sample(max_signatures, random_state=seed)
    genes = read_gene_info(gene_info)
    lm = genes[genes["pr_is_lm"].eq(1)]["gene_id"].astype(str).tolist()
    all_ids = genes["gene_id"].astype(str).tolist()

    x_lm = read_gctx_block(gctx, sig["sig_id"].astype(str).tolist(), lm, chunk_size=chunk_size)
    x_mean = x_lm.mean(axis=0)
    x_center = x_lm - x_mean
    xtx = x_center.T @ x_center
    xtx.flat[:: xtx.shape[0] + 1] += alpha
    inv_xtx = np.linalg.inv(xtx).astype(np.float32)

    coefs = np.empty((len(lm), len(all_ids)), dtype=np.float32)
    y_mean = np.empty(len(all_ids), dtype=np.float32)
    for start in range(0, len(all_ids), 512):
        stop = min(start + 512, len(all_ids))
        y = read_gctx_block(gctx, sig["sig_id"].astype(str).tolist(), all_ids[start:stop], chunk_size=chunk_size)
        ym = y.mean(axis=0)
        y_mean[start:stop] = ym
        coefs[:, start:stop] = inv_xtx @ (x_center.T @ (y - ym))

    out_npz = Path(out_npz)
    out_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out_npz,
        coef=coefs,
        x_mean=x_mean.astype(np.float32),
        y_mean=y_mean.astype(np.float32),
        landmark_gene_ids=np.array(lm),
        all_gene_ids=np.array(all_ids),
        all_gene_names=genes["gene_name"].astype(str).to_numpy(),
        alpha=np.array([alpha], dtype=np.float32),
    )


def apply_landmark_model(
    landmark_pred_h5ad: str | Path,
    model_npz: str | Path,
    out_h5ad: str | Path,
) -> None:
    import anndata as ad
    import scanpy as sc

    pred = sc.read_h5ad(landmark_pred_h5ad)
    model = np.load(model_npz, allow_pickle=True)
    landmark_ids = model["landmark_gene_ids"].astype(str)
    all_names = model["all_gene_names"].astype(str)
    pred_gene_ids = pred.var["gene_id"].astype(str) if "gene_id" in pred.var else pred.var_names.astype(str)
    keep = [i for i, gid in enumerate(landmark_ids) if gid in set(pred_gene_ids)]
    if len(keep) != len(landmark_ids):
        missing = len(landmark_ids) - len(keep)
        raise ValueError(f"Landmark prediction matrix is missing {missing} required landmark genes")
    order = [np.where(pred_gene_ids == gid)[0][0] for gid in landmark_ids]
    x = np.asarray(pred.X[:, order], dtype=np.float32)
    y = (x - model["x_mean"]) @ model["coef"] + model["y_mean"]
    var = pd.DataFrame(
        {
            "gene_id": model["all_gene_ids"].astype(str),
            "gene_name": all_names,
            "gene_response_source": np.where(
                np.isin(model["all_gene_ids"].astype(str), landmark_ids),
                "landmark_measured",
                "post_cpa_imputed",
            ),
        },
        index=all_names,
    )
    out = ad.AnnData(X=y.astype(np.float32), obs=pred.obs.copy(), var=var)
    out_h5ad = Path(out_h5ad)
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    out.write_h5ad(out_h5ad, compression="gzip")

