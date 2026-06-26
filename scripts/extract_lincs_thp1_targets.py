#!/usr/bin/env python3
"""Extract THP1 compound signatures for coloc target genes from LINCS GCTX."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm


RAW = Path("data/raw/lincs_gse92742")
DEFAULT_GCTX = RAW / "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
DEFAULT_SIG = RAW / "GSE92742_Broad_LINCS_sig_info.txt.gz"
DEFAULT_GENE = RAW / "GSE92742_Broad_LINCS_gene_info.txt.gz"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", default="data/processed/protective_expression_gene_summary.tsv")
    parser.add_argument("--gctx", default=str(DEFAULT_GCTX))
    parser.add_argument("--sig-info", default=str(DEFAULT_SIG))
    parser.add_argument("--gene-info", default=str(DEFAULT_GENE))
    parser.add_argument("--cell-id", default="THP1")
    parser.add_argument("--pert-type", default="trt_cp")
    parser.add_argument("--out-long", default="data/processed/lincs_thp1_target_gene_zscores.tsv.gz")
    parser.add_argument("--out-summary", default="data/processed/lincs_thp1_drug_gene_summary.tsv")
    parser.add_argument("--out-coverage", default="data/processed/lincs_target_gene_coverage.tsv")
    return parser.parse_args()


def decode_array(arr: np.ndarray) -> list[str]:
    out = []
    for x in arr:
        if isinstance(x, bytes):
            out.append(x.decode("utf-8"))
        else:
            out.append(str(x))
    return out


def gctx_handles(h5: h5py.File) -> tuple[h5py.Dataset, list[str], list[str], str]:
    matrix = h5["0/DATA/0/matrix"]
    gene_ids = decode_array(h5["0/META/ROW/id"][:])
    sig_ids = decode_array(h5["0/META/COL/id"][:])
    if matrix.shape == (len(gene_ids), len(sig_ids)):
        return matrix, gene_ids, sig_ids, "row"
    if matrix.shape == (len(sig_ids), len(gene_ids)):
        return matrix, gene_ids, sig_ids, "col"
    raise ValueError(
        f"Could not reconcile GCTX matrix shape {matrix.shape} with "
        f"{len(gene_ids)} row IDs and {len(sig_ids)} column IDs"
    )


def read_gene_info(path: str) -> pd.DataFrame:
    gene = pd.read_csv(path, sep="\t", compression="infer", dtype=str)
    symbol_col = "pr_gene_symbol" if "pr_gene_symbol" in gene.columns else "gene_symbol"
    id_col = "pr_gene_id" if "pr_gene_id" in gene.columns else gene.columns[0]
    return gene.rename(columns={symbol_col: "gene_name", id_col: "gene_id"})


def read_sig_info(path: str, cell_id: str, pert_type: str) -> pd.DataFrame:
    sig = pd.read_csv(path, sep="\t", compression="infer", dtype=str)
    keep = sig["cell_id"].eq(cell_id) & sig["pert_type"].eq(pert_type)
    sig = sig.loc[keep].copy()
    rename = {}
    if "pert_idose" in sig.columns and "pert_dose" not in sig.columns:
        rename["pert_idose"] = "pert_dose"
    if "pert_itime" in sig.columns and "pert_time" not in sig.columns:
        rename["pert_itime"] = "pert_time"
    if rename:
        sig = sig.rename(columns=rename)
    for col in ["pert_dose", "pert_time", "distil_cc_q75", "pct_self_rank_q25"]:
        if col in sig.columns:
            sig[col] = pd.to_numeric(sig[col].astype(str).str.extract(r"([-+]?[0-9]*\.?[0-9]+)", expand=False), errors="coerce")
    return sig


def main() -> None:
    args = parse_args()
    targets = pd.read_csv(args.targets, sep="\t")
    target_genes = sorted(targets["gene_name"].dropna().astype(str).unique())
    gene_info = read_gene_info(args.gene_info)
    target_gene_info = gene_info.loc[gene_info["gene_name"].isin(target_genes), ["gene_id", "gene_name"]].drop_duplicates()
    target_gene_info["in_lincs"] = True

    coverage = pd.DataFrame({"gene_name": target_genes}).merge(target_gene_info, on="gene_name", how="left")
    coverage["in_lincs"] = coverage["in_lincs"].fillna(False)
    Path(args.out_coverage).parent.mkdir(parents=True, exist_ok=True)
    coverage.to_csv(args.out_coverage, sep="\t", index=False)

    sig = read_sig_info(args.sig_info, args.cell_id, args.pert_type)
    if sig.empty:
        raise ValueError(f"No signatures matched cell_id={args.cell_id!r}, pert_type={args.pert_type!r}")

    with h5py.File(args.gctx, "r") as h5:
        matrix, gene_ids, sig_ids, gene_axis = gctx_handles(h5)
        gene_lookup = {rid: i for i, rid in enumerate(gene_ids)}
        sig_lookup = {cid: i for i, cid in enumerate(sig_ids)}

        row_meta = target_gene_info[target_gene_info["gene_id"].astype(str).isin(gene_lookup)].copy()
        sig = sig[sig["sig_id"].isin(sig_lookup)].copy()
        if row_meta.empty:
            raise ValueError("None of the target genes were found in the GCTX row IDs")
        if sig.empty:
            raise ValueError("None of the THP1 signature IDs were found in the GCTX column IDs")

        sig_order = np.array([sig_lookup[x] for x in sig["sig_id"]])
        sig_sort = np.argsort(sig_order)
        sorted_sigs = sig_order[sig_sort]
        sorted_sig = sig.iloc[sig_sort].reset_index(drop=True)

        records = []
        for _, gene_row in tqdm(row_meta.iterrows(), total=len(row_meta), desc="Extracting genes"):
            rid = str(gene_row["gene_id"])
            gene_idx = gene_lookup[rid]
            if gene_axis == "row":
                values = matrix[gene_idx, sorted_sigs]
            else:
                values = matrix[sorted_sigs, gene_idx]
            frame = sorted_sig.copy()
            frame["gene_id"] = rid
            frame["gene_name"] = gene_row["gene_name"]
            frame["zscore"] = np.asarray(values, dtype=float)
            records.append(frame)

    long = pd.concat(records, ignore_index=True)
    meta_cols = [
        c
        for c in [
            "sig_id",
            "pert_id",
            "pert_iname",
            "pert_type",
            "cell_id",
            "pert_dose",
            "pert_dose_unit",
            "pert_time",
            "pert_time_unit",
            "distil_cc_q75",
            "pct_self_rank_q25",
            "gene_id",
            "gene_name",
            "zscore",
        ]
        if c in long.columns
    ]
    long = long[meta_cols]
    with gzip.open(args.out_long, "wt") as fout:
        long.to_csv(fout, sep="\t", index=False)

    summary = (
        long.groupby(["pert_id", "pert_iname", "gene_name"], dropna=False, as_index=False)
        .agg(mean_z=("zscore", "mean"), median_z=("zscore", "median"), n_signatures=("sig_id", "nunique"))
        .sort_values(["pert_iname", "gene_name"])
    )
    summary.to_csv(args.out_summary, sep="\t", index=False)

    print(f"Wrote {len(long):,} gene-signature rows to {args.out_long}")
    print(f"Wrote {len(summary):,} drug-gene rows to {args.out_summary}")
    print(f"Wrote LINCS target coverage to {args.out_coverage}")


if __name__ == "__main__":
    main()
