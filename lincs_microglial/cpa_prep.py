"""Prepare LINCS Level 5 data for CPA."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .lincs_io import (
    DEFAULT_GCTX,
    DEFAULT_GENE_INFO,
    DEFAULT_PERT_INFO,
    DEFAULT_SIG_INFO,
    read_gene_info,
    read_gctx_block,
    read_pert_info,
    read_sig_info,
    select_gene_info,
)


CONTROL_GROUP = "DMSO"
CPA_DOSE_KEY = "cpa_dose"


@dataclass(frozen=True)
class PreparedCpaData:
    obs: pd.DataFrame
    var: pd.DataFrame
    x: np.ndarray
    eligible_drugs: pd.DataFrame


def read_target_genes(path: str | Path | None) -> list[str]:
    if path is None or not Path(path).exists():
        return []
    targets = pd.read_csv(path, sep="\t")
    return sorted(targets["gene_name"].dropna().astype(str).unique())


def add_drug_metadata(sig: pd.DataFrame, pert: pd.DataFrame) -> pd.DataFrame:
    cols = ["pert_id", "canonical_smiles", "pubchem_cid", "inchi_key", "inchi_key_prefix"]
    cols = [c for c in cols if c in pert.columns]
    out = sig.merge(pert[cols].drop_duplicates("pert_id"), on="pert_id", how="left")
    out["canonical_smiles"] = out["canonical_smiles"].replace({"-666": np.nan, "": np.nan})
    return out


def eligible_drug_table(
    sig: pd.DataFrame,
    min_non_thp1_signatures: int = 20,
    min_non_thp1_cell_lines: int = 3,
    thp1_cell_id: str = "THP1",
) -> pd.DataFrame:
    trt = sig[sig["pert_type"].eq("trt_cp")].copy()
    non = trt[~trt["cell_id"].eq(thp1_cell_id)]
    grouped = (
        non.groupby(["pert_id", "pert_iname"], as_index=False)
        .agg(
            n_non_thp1_signatures=("sig_id", "nunique"),
            n_non_thp1_cell_lines=("cell_id", "nunique"),
            n_non_thp1_doses=("pert_dose", "nunique"),
            canonical_smiles=("canonical_smiles", first_nonmissing),
        )
    )
    thp1 = (
        trt[trt["cell_id"].eq(thp1_cell_id)]
        .groupby("pert_id", as_index=False)
        .agg(n_thp1_signatures=("sig_id", "nunique"), n_thp1_doses=("pert_dose", "nunique"))
    )
    out = grouped.merge(thp1, on="pert_id", how="left")
    out[["n_thp1_signatures", "n_thp1_doses"]] = out[["n_thp1_signatures", "n_thp1_doses"]].fillna(0).astype(int)
    out["has_thp1_observed"] = out["n_thp1_signatures"].gt(0)
    out["has_smiles"] = out["canonical_smiles"].notna()
    out["eligible"] = (
        out["has_smiles"]
        & out["n_non_thp1_signatures"].ge(min_non_thp1_signatures)
        & out["n_non_thp1_cell_lines"].ge(min_non_thp1_cell_lines)
    )
    return out.sort_values(["eligible", "n_non_thp1_signatures"], ascending=[False, False])


def first_nonmissing(values: pd.Series) -> str | float:
    vals = values.dropna().astype(str)
    return vals.iloc[0] if len(vals) else np.nan


def add_cpa_columns(sig: pd.DataFrame) -> pd.DataFrame:
    out = sig.copy()
    out["condition_ID"] = np.where(out["pert_type"].str.startswith("ctl_"), CONTROL_GROUP, out["pert_id"])
    out["smiles_rdkit"] = np.where(out["condition_ID"].eq(CONTROL_GROUP), "", out["canonical_smiles"].fillna(""))
    dose = pd.to_numeric(out["pert_dose"], errors="coerce")
    out["dose_um"] = dose
    out["log_dose"] = np.where(out["condition_ID"].eq(CONTROL_GROUP), 0.0, np.log10(dose.clip(lower=1e-6)))
    out["log_dose"] = pd.Series(out["log_dose"]).replace([np.inf, -np.inf], np.nan).fillna(0.0).astype(float)
    cpa_dose = np.where(out["condition_ID"].eq(CONTROL_GROUP), 0.0, dose.clip(lower=1e-6))
    out[CPA_DOSE_KEY] = pd.Series(cpa_dose).replace([np.inf, -np.inf], np.nan).fillna(0.0).astype(float)
    out["cell_type"] = out["cell_id"].astype(str)
    return out


def assign_splits(
    obs: pd.DataFrame,
    mode: str,
    thp1_cell_id: str = "THP1",
    valid_fraction: float = 0.1,
    seed: int = 1,
) -> pd.Series:
    rng = np.random.default_rng(seed)
    split = pd.Series("train", index=obs.index, dtype=object)
    is_query = obs["row_kind"].eq("query")
    is_thp1_trt = obs["cell_id"].eq(thp1_cell_id) & obs["pert_type"].eq("trt_cp") & ~is_query
    split.loc[is_query] = "query"
    if mode == "validation":
        split.loc[is_thp1_trt] = "ood"
    elif mode == "final":
        split.loc[is_thp1_trt] = "train"
    else:
        raise ValueError("mode must be 'validation' or 'final'")

    trainable = split.eq("train") & obs["pert_type"].eq("trt_cp") & ~obs["cell_id"].eq(thp1_cell_id)
    combos = obs.loc[trainable, ["pert_id", "cell_id"]].drop_duplicates()
    if len(combos) and valid_fraction > 0:
        n_valid = max(1, int(round(len(combos) * valid_fraction)))
        valid_idx = rng.choice(combos.index.to_numpy(), size=min(n_valid, len(combos)), replace=False)
        valid_combos = set(map(tuple, combos.loc[valid_idx, ["pert_id", "cell_id"]].to_numpy()))
        valid_rows = trainable & obs[["pert_id", "cell_id"]].apply(tuple, axis=1).isin(valid_combos)
        split.loc[valid_rows] = "test"
    return split


def query_rows(
    obs: pd.DataFrame,
    eligible: pd.DataFrame,
    mode: str,
    thp1_cell_id: str = "THP1",
) -> pd.DataFrame:
    controls = obs[obs["cell_id"].eq(thp1_cell_id) & obs["condition_ID"].eq(CONTROL_GROUP)].copy()
    if controls.empty:
        raise ValueError(f"No {thp1_cell_id} vehicle controls are available for CPA query rows")
    base = controls.iloc[[0]].copy()

    eligible = eligible[eligible["eligible"]].copy()
    if mode == "final":
        eligible = eligible[~eligible["has_thp1_observed"]].copy()
    rows = []
    dose_lookup = (
        obs[obs["pert_type"].eq("trt_cp")]
        .dropna(subset=["pert_dose"])
        .groupby("pert_id")["pert_dose"]
        .agg(lambda x: float(pd.Series(x).mode().iloc[0]) if len(pd.Series(x).mode()) else float(np.nanmedian(x)))
    )
    for drug in eligible.itertuples(index=False):
        row = base.copy()
        row["sig_id"] = f"CPA_QUERY_THP1:{drug.pert_id}"
        row["pert_id"] = drug.pert_id
        row["pert_iname"] = drug.pert_iname
        row["pert_type"] = "trt_cp"
        row["cell_id"] = thp1_cell_id
        row["cell_type"] = thp1_cell_id
        row["condition_ID"] = drug.pert_id
        row["canonical_smiles"] = drug.canonical_smiles
        row["smiles_rdkit"] = drug.canonical_smiles
        dose = dose_lookup.get(drug.pert_id, 10.0)
        row["pert_dose"] = dose
        row["dose_um"] = dose
        row["log_dose"] = float(np.log10(max(dose, 1e-6)))
        row[CPA_DOSE_KEY] = float(max(dose, 1e-6))
        row["row_kind"] = "query"
        rows.append(row)
    if not rows:
        return pd.DataFrame(columns=obs.columns)
    return pd.concat(rows, ignore_index=True)


def prepare_cpa_data(
    gctx: str | Path = DEFAULT_GCTX,
    sig_info: str | Path = DEFAULT_SIG_INFO,
    gene_info: str | Path = DEFAULT_GENE_INFO,
    pert_info: str | Path = DEFAULT_PERT_INFO,
    targets: str | Path | None = "data/processed/protective_expression_gene_summary.tsv",
    mode: str = "validation",
    gene_mode: str = "pilot",
    time_hours: float = 6,
    min_non_thp1_signatures: int = 20,
    min_non_thp1_cell_lines: int = 3,
    thp1_cell_id: str = "THP1",
    seed: int = 1,
    chunk_size: int = 2048,
) -> PreparedCpaData:
    sig = read_sig_info(sig_info)
    pert = read_pert_info(pert_info)
    sig = add_drug_metadata(sig, pert)
    sig = sig[sig["pert_time"].eq(time_hours)].copy()
    is_control = sig["pert_type"].str.startswith("ctl_vehicle")
    is_trt = sig["pert_type"].eq("trt_cp")
    sig = sig[is_control | is_trt].copy()

    eligible = eligible_drug_table(sig, min_non_thp1_signatures, min_non_thp1_cell_lines, thp1_cell_id)
    eligible_ids = set(eligible.loc[eligible["eligible"], "pert_id"])
    keep = sig["pert_type"].str.startswith("ctl_vehicle") | sig["pert_id"].isin(eligible_ids)
    obs = sig.loc[keep].copy()
    obs["row_kind"] = "observed"
    obs = add_cpa_columns(obs)
    q = query_rows(obs, eligible, mode=mode, thp1_cell_id=thp1_cell_id)
    obs = pd.concat([obs, q], ignore_index=True)
    obs["split"] = assign_splits(obs, mode=mode, thp1_cell_id=thp1_cell_id, seed=seed)

    observed_obs = obs[obs["row_kind"].eq("observed")].copy()
    target_genes = read_target_genes(targets)
    var = select_gene_info(read_gene_info(gene_info), gene_mode, target_genes=target_genes, n_variable=2000)
    x_obs = read_gctx_block(gctx, observed_obs["sig_id"].astype(str).tolist(), var["gene_id"].astype(str).tolist(), chunk_size)
    if q.empty:
        x = x_obs
    else:
        thp1_controls = np.where(
            observed_obs["cell_id"].eq(thp1_cell_id).to_numpy()
            & observed_obs["condition_ID"].eq(CONTROL_GROUP).to_numpy()
        )[0]
        if len(thp1_controls) == 0:
            raise ValueError("Cannot build query rows without observed THP1 control rows")
        control_mean = x_obs[thp1_controls, :].mean(axis=0, keepdims=True).astype(np.float32)
        x_query = np.repeat(control_mean, repeats=len(q), axis=0)
        x = np.vstack([x_obs, x_query]).astype(np.float32)
    obs = obs.reset_index(drop=True)
    obs.index = obs["sig_id"].astype(str)
    var = var.reset_index(drop=True)
    var.index = var["gene_name"].astype(str)
    return PreparedCpaData(obs=obs, var=var, x=x, eligible_drugs=eligible)


def write_anndata(prepared: PreparedCpaData, out_h5ad: str | Path, out_eligible: str | Path | None = None) -> None:
    try:
        import anndata as ad
    except ImportError as exc:
        raise ImportError("anndata is required to write CPA input; install/use the CPA environment") from exc

    out_h5ad = Path(out_h5ad)
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata = ad.AnnData(X=prepared.x, obs=prepared.obs, var=prepared.var)
    adata.write_h5ad(out_h5ad, compression="gzip")
    if out_eligible is not None:
        out_eligible = Path(out_eligible)
        out_eligible.parent.mkdir(parents=True, exist_ok=True)
        prepared.eligible_drugs.to_csv(out_eligible, sep="\t", index=False)
