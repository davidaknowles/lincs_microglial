#!/usr/bin/env python3
"""Summarize CPA validation predictions against observed THP1 signatures."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observed-gene", default="data/processed/lincs_thp1_drug_gene_summary.tsv")
    parser.add_argument("--predicted-gene", required=True)
    parser.add_argument("--out", default="data/processed/cpa/cpa_validation_metrics.tsv")
    return parser.parse_args()


def corr(x: pd.Series, y: pd.Series) -> float:
    if len(x) < 2 or x.nunique() < 2 or y.nunique() < 2:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def main() -> None:
    args = parse_args()
    obs = pd.read_csv(args.observed_gene, sep="\t")
    pred = pd.read_csv(args.predicted_gene, sep="\t")
    merged = pred.merge(
        obs[["pert_id", "gene_name", "mean_z"]].rename(columns={"mean_z": "observed_mean_z"}),
        on=["pert_id", "gene_name"],
        how="inner",
    )
    merged = merged.rename(columns={"mean_z": "predicted_mean_z"})
    if merged.empty:
        raise ValueError("No overlapping observed/predicted drug-gene rows for validation")
    rows = []
    for (pert_id, pert_iname), frame in merged.groupby(["pert_id", "pert_iname"], dropna=False):
        rows.append(
            {
                "pert_id": pert_id,
                "pert_iname": pert_iname,
                "n_genes": frame["gene_name"].nunique(),
                "pearson_r": corr(frame["predicted_mean_z"], frame["observed_mean_z"]),
                "mae": float(np.mean(np.abs(frame["predicted_mean_z"] - frame["observed_mean_z"]))),
                "rmse": float(np.sqrt(np.mean((frame["predicted_mean_z"] - frame["observed_mean_z"]) ** 2))),
            }
        )
    by_drug = pd.DataFrame(rows)
    overall = pd.DataFrame(
        {
            "pert_id": ["OVERALL"],
            "pert_iname": ["OVERALL"],
            "n_genes": [merged["gene_name"].nunique()],
            "pearson_r": [corr(merged["predicted_mean_z"], merged["observed_mean_z"])],
            "mae": [float(np.mean(np.abs(merged["predicted_mean_z"] - merged["observed_mean_z"])))],
            "rmse": [float(np.sqrt(np.mean((merged["predicted_mean_z"] - merged["observed_mean_z"]) ** 2)))],
        }
    )
    out = pd.concat([overall, by_drug], ignore_index=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote CPA validation metrics: {args.out}")


if __name__ == "__main__":
    main()
