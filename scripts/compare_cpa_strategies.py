#!/usr/bin/env python3
"""Compare two CPA THP1 prediction strategies on common drugs and genes."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observed-gene", default="data/processed/lincs_thp1_drug_gene_summary.tsv")
    parser.add_argument("--a-gene", required=True)
    parser.add_argument("--b-gene", required=True)
    parser.add_argument("--a-name", default="strategy_a")
    parser.add_argument("--b-name", default="strategy_b")
    parser.add_argument("--a-drug-scores", default=None)
    parser.add_argument("--b-drug-scores", default=None)
    parser.add_argument("--out-metrics", required=True)
    parser.add_argument("--out-common-gene", required=True)
    parser.add_argument("--out-rank", default=None)
    return parser.parse_args()


def corr(x: pd.Series, y: pd.Series) -> float:
    if len(x) < 2 or x.nunique() < 2 or y.nunique() < 2:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def metrics(label: str, frame: pd.DataFrame) -> dict[str, object]:
    diff = frame[f"{label}_mean_z"] - frame["observed_mean_z"]
    return {
        "strategy": label,
        "n_drugs": frame["pert_id"].nunique(),
        "n_genes": frame["gene_name"].nunique(),
        "n_drug_gene_rows": len(frame),
        "pearson_r": corr(frame[f"{label}_mean_z"], frame["observed_mean_z"]),
        "mae": float(np.mean(np.abs(diff))),
        "rmse": float(np.sqrt(np.mean(diff**2))),
    }


def read_gene(path: str, label: str) -> pd.DataFrame:
    cols = ["pert_id", "pert_iname", "gene_name", "mean_z"]
    return pd.read_csv(path, sep="\t", usecols=lambda c: c in cols).rename(columns={"mean_z": f"{label}_mean_z"})


def main() -> None:
    args = parse_args()
    obs = pd.read_csv(args.observed_gene, sep="\t", usecols=lambda c: c in {"pert_id", "gene_name", "mean_z"})
    obs = obs.rename(columns={"mean_z": "observed_mean_z"})
    a = read_gene(args.a_gene, args.a_name)
    b = read_gene(args.b_gene, args.b_name)
    common = (
        a.merge(b[["pert_id", "gene_name", f"{args.b_name}_mean_z"]], on=["pert_id", "gene_name"], how="inner")
        .merge(obs, on=["pert_id", "gene_name"], how="inner")
        .sort_values(["pert_iname", "gene_name"])
    )
    if common.empty:
        raise ValueError("No common observed/predicted drug-gene rows")

    rows = [metrics(args.a_name, common), metrics(args.b_name, common)]
    out_metrics = pd.DataFrame(rows)
    Path(args.out_metrics).parent.mkdir(parents=True, exist_ok=True)
    out_metrics.to_csv(args.out_metrics, sep="\t", index=False)
    Path(args.out_common_gene).parent.mkdir(parents=True, exist_ok=True)
    common.to_csv(args.out_common_gene, sep="\t", index=False)

    if args.a_drug_scores and args.b_drug_scores and args.out_rank:
        score_cols = ["pert_id", "pert_iname", "n_protective_genes", "n_opposing_genes", "frac_protective_genes"]
        da = pd.read_csv(args.a_drug_scores, sep="\t", usecols=lambda c: c in score_cols).copy()
        db = pd.read_csv(args.b_drug_scores, sep="\t", usecols=lambda c: c in score_cols).copy()
        da["a_rank"] = np.arange(1, len(da) + 1)
        db["b_rank"] = np.arange(1, len(db) + 1)
        rank = da.merge(db, on=["pert_id", "pert_iname"], suffixes=(f"_{args.a_name}", f"_{args.b_name}"))
        rank["rank_delta_b_minus_a"] = rank["b_rank"] - rank["a_rank"]
        Path(args.out_rank).parent.mkdir(parents=True, exist_ok=True)
        rank.to_csv(args.out_rank, sep="\t", index=False)

    print(f"Wrote strategy metrics: {args.out_metrics}")
    print(f"Wrote common drug-gene table: {args.out_common_gene}")


if __name__ == "__main__":
    main()
