"""LINCS GSE92742 metadata and GCTX helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


RAW = Path("data/raw/lincs_gse92742")
DEFAULT_GCTX = RAW / "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
DEFAULT_SIG_INFO = RAW / "GSE92742_Broad_LINCS_sig_info.txt.gz"
DEFAULT_GENE_INFO = RAW / "GSE92742_Broad_LINCS_gene_info.txt.gz"
DEFAULT_PERT_INFO = RAW / "GSE92742_Broad_LINCS_pert_info.txt.gz"


@dataclass(frozen=True)
class GctxIndex:
    """Resolved GCTX matrix axes."""

    matrix: h5py.Dataset
    gene_ids: list[str]
    sig_ids: list[str]
    gene_axis: str


def decode_array(arr: np.ndarray) -> list[str]:
    out: list[str] = []
    for x in arr:
        out.append(x.decode("utf-8") if isinstance(x, bytes) else str(x))
    return out


def gctx_index(h5: h5py.File) -> GctxIndex:
    matrix = h5["0/DATA/0/matrix"]
    gene_ids = decode_array(h5["0/META/ROW/id"][:])
    sig_ids = decode_array(h5["0/META/COL/id"][:])
    if matrix.shape == (len(gene_ids), len(sig_ids)):
        return GctxIndex(matrix=matrix, gene_ids=gene_ids, sig_ids=sig_ids, gene_axis="row")
    if matrix.shape == (len(sig_ids), len(gene_ids)):
        return GctxIndex(matrix=matrix, gene_ids=gene_ids, sig_ids=sig_ids, gene_axis="col")
    raise ValueError(
        f"Could not reconcile GCTX matrix shape {matrix.shape} with "
        f"{len(gene_ids)} gene IDs and {len(sig_ids)} signature IDs"
    )


def numeric_from_lincs(value: pd.Series) -> pd.Series:
    """Parse LINCS dose/time columns that may contain values like '10 uM'."""

    return pd.to_numeric(
        value.astype(str).str.extract(r"([-+]?[0-9]*\.?[0-9]+)", expand=False),
        errors="coerce",
    )


def read_gene_info(path: str | Path = DEFAULT_GENE_INFO) -> pd.DataFrame:
    gene = pd.read_csv(path, sep="\t", compression="infer", dtype=str)
    rename = {}
    if "pr_gene_id" in gene.columns:
        rename["pr_gene_id"] = "gene_id"
    if "pr_gene_symbol" in gene.columns:
        rename["pr_gene_symbol"] = "gene_name"
    gene = gene.rename(columns=rename)
    gene["gene_id"] = gene["gene_id"].astype(str)
    gene["gene_name"] = gene["gene_name"].astype(str)
    for col in ["pr_is_lm", "pr_is_bing"]:
        if col in gene.columns:
            gene[col] = pd.to_numeric(gene[col], errors="coerce").fillna(0).astype(int)
    return gene


def read_sig_info(path: str | Path = DEFAULT_SIG_INFO) -> pd.DataFrame:
    sig = pd.read_csv(path, sep="\t", compression="infer", dtype=str)
    for col in ["pert_dose", "pert_time"]:
        if col in sig.columns:
            sig[col] = numeric_from_lincs(sig[col])
    return sig


def read_pert_info(path: str | Path = DEFAULT_PERT_INFO) -> pd.DataFrame:
    pert = pd.read_csv(path, sep="\t", compression="infer", dtype=str)
    for col in ["canonical_smiles", "pert_iname"]:
        if col in pert.columns:
            pert[col] = pert[col].replace({"-666": np.nan, "": np.nan})
    return pert


def select_gene_info(
    gene_info: pd.DataFrame,
    mode: str,
    target_genes: list[str] | None = None,
    n_variable: int | None = None,
) -> pd.DataFrame:
    """Select genes for CPA input.

    `n_variable` is a deterministic placeholder for pilot-sized runs when no
    variance table exists yet; it keeps landmarks first, then fills by GCTX order.
    """

    gene = gene_info.copy()
    target_genes = set(target_genes or [])
    if mode == "all":
        keep = pd.Series(True, index=gene.index)
    elif mode == "landmark":
        keep = gene.get("pr_is_lm", 0).eq(1)
    elif mode == "pilot":
        if n_variable is None:
            n_variable = 2000
        landmark = gene.get("pr_is_lm", 0).eq(1)
        keep = landmark | gene["gene_name"].isin(target_genes)
        if keep.sum() < n_variable:
            fill = gene.index[~keep][: max(0, n_variable - int(keep.sum()))]
            keep.loc[fill] = True
    else:
        raise ValueError(f"Unknown gene mode: {mode}")
    if target_genes:
        keep |= gene["gene_name"].isin(target_genes)
    return gene.loc[keep, ["gene_id", "gene_name", "pr_is_lm", "pr_is_bing"]].drop_duplicates("gene_id")


def read_gctx_block(
    gctx_path: str | Path,
    sig_ids: list[str],
    gene_ids: list[str],
    chunk_size: int = 2048,
) -> np.ndarray:
    """Read a signature x gene dense block from a GCTX file."""

    with h5py.File(gctx_path, "r") as h5:
        idx = gctx_index(h5)
        sig_lookup = {sid: i for i, sid in enumerate(idx.sig_ids)}
        gene_lookup = {gid: i for i, gid in enumerate(idx.gene_ids)}
        missing_sigs = [x for x in sig_ids if x not in sig_lookup]
        missing_genes = [x for x in gene_ids if x not in gene_lookup]
        if missing_sigs:
            raise ValueError(f"{len(missing_sigs)} requested signatures are absent from GCTX")
        if missing_genes:
            raise ValueError(f"{len(missing_genes)} requested genes are absent from GCTX")

        sig_pos = np.array([sig_lookup[x] for x in sig_ids])
        gene_pos = np.array([gene_lookup[x] for x in gene_ids])
        sig_sort = np.argsort(sig_pos)
        sorted_sig_pos = sig_pos[sig_sort]
        reverse_sig = np.argsort(sig_sort)
        x_sorted = np.empty((len(sig_ids), len(gene_ids)), dtype=np.float32)

        for start in range(0, len(sorted_sig_pos), chunk_size):
            stop = min(start + chunk_size, len(sorted_sig_pos))
            sig_chunk = sorted_sig_pos[start:stop]
            if idx.gene_axis == "row":
                block = np.asarray(idx.matrix[:, sig_chunk], dtype=np.float32).T
                block = block[:, gene_pos]
            else:
                block = np.asarray(idx.matrix[sig_chunk, :], dtype=np.float32)
                block = block[:, gene_pos]
            x_sorted[start:stop, :] = block
        return x_sorted[reverse_sig, :]

