#!/usr/bin/env python3
"""Score THP1 LINCS compounds by protective push across ISOMIGA target genes."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", default="data/processed/protective_expression_gene_summary.tsv")
    parser.add_argument("--drug-gene", default="data/processed/lincs_thp1_drug_gene_summary.tsv")
    parser.add_argument("--out-drug", default="data/processed/lincs_thp1_protective_drug_scores.tsv")
    parser.add_argument("--out-gene", default="data/processed/lincs_thp1_protective_drug_gene_scores.tsv")
    parser.add_argument("--min-genes", type=int, default=2)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    targets = pd.read_csv(args.targets, sep="\t")
    drug_gene = pd.read_csv(args.drug_gene, sep="\t")

    target_cols = [
        "gene_name",
        "protective_expression_direction",
        "protective_direction_label",
        "max_h4",
        "mean_h4",
        "n_coloc",
        "total_weight",
        "max_gwas_z_abs",
        "mr_ivw_beta",
        "mr_ivw_se",
        "mr_ivw_p",
    ]
    merged = drug_gene.merge(targets[target_cols], on="gene_name", how="inner")
    merged["protective_push_z"] = merged["mean_z"] * merged["protective_expression_direction"]
    merged["weighted_protective_push"] = merged["protective_push_z"] * merged["total_weight"]
    merged["pushes_protective"] = merged["protective_push_z"].gt(0)
    merged["strong_push"] = merged["protective_push_z"].ge(1.0)
    merged["opposes_protective"] = merged["protective_push_z"].lt(0)

    denom = merged.groupby(["pert_id", "pert_iname"])["total_weight"].transform("sum")
    merged["gene_weight_fraction"] = np.where(denom > 0, merged["total_weight"] / denom, np.nan)

    drug = (
        merged.groupby(["pert_id", "pert_iname"], as_index=False)
        .agg(
            n_target_genes=("gene_name", "nunique"),
            n_protective_genes=("pushes_protective", "sum"),
            n_strong_protective_genes=("strong_push", "sum"),
            n_opposing_genes=("opposes_protective", "sum"),
            mean_protective_push_z=("protective_push_z", "mean"),
            median_protective_push_z=("protective_push_z", "median"),
            summed_weighted_push=("weighted_protective_push", "sum"),
            total_weight=("total_weight", "sum"),
            mean_target_mr_ivw_beta=("mr_ivw_beta", "mean"),
            median_target_mr_ivw_beta=("mr_ivw_beta", "median"),
            min_target_mr_ivw_p=("mr_ivw_p", "min"),
            min_n_signatures=("n_signatures", "min"),
            total_gene_signatures=("n_signatures", "sum"),
        )
    )
    drug["weighted_mean_protective_push_z"] = drug["summed_weighted_push"] / drug["total_weight"].replace(0, np.nan)
    drug["fraction_genes_protective"] = drug["n_protective_genes"] / drug["n_target_genes"].replace(0, np.nan)
    drug = drug[drug["n_target_genes"].ge(args.min_genes)].copy()
    drug = drug.sort_values(
        ["n_strong_protective_genes", "weighted_mean_protective_push_z", "fraction_genes_protective", "n_target_genes"],
        ascending=False,
    )

    Path(args.out_drug).parent.mkdir(parents=True, exist_ok=True)
    drug.to_csv(args.out_drug, sep="\t", index=False)
    merged.sort_values(["pert_iname", "protective_push_z"], ascending=[True, False]).to_csv(
        args.out_gene, sep="\t", index=False
    )

    print(f"Wrote {len(drug):,} drug scores to {args.out_drug}")
    print(f"Wrote {len(merged):,} drug-gene scores to {args.out_gene}")


if __name__ == "__main__":
    main()
