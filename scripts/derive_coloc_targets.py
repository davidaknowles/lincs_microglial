#!/usr/bin/env python3
"""Derive protective expression targets from ISOMIGA AD colocalization results."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coloc", default="isomiga_AD_coloc.txt", help="ISOMIGA coloc TSV")
    parser.add_argument("--out", default="data/processed/protective_expression_targets.tsv")
    parser.add_argument("--summary-out", default="data/processed/protective_expression_gene_summary.tsv")
    parser.add_argument("--h4-threshold", type=float, default=0.8)
    parser.add_argument("--require-pass-distance", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument(
        "--set-label",
        default="expression",
        help="Case-insensitive substring to keep in set_label. Use empty string to disable.",
    )
    return parser.parse_args()


def sign_label(x: float) -> str:
    if x > 0:
        return "increase"
    if x < 0:
        return "decrease"
    return "neutral"


def main() -> None:
    args = parse_args()
    coloc = pd.read_csv(args.coloc, sep="\t")

    numeric_cols = ["QTL_Beta", "QTL_SE", "GWAS_SNP_Beta", "GWAS_SNP_SE", "PP.H4.abf", "GWAS_P"]
    for col in numeric_cols:
        coloc[col] = pd.to_numeric(coloc[col], errors="coerce")

    keep = coloc["PP.H4.abf"].ge(args.h4_threshold)
    if args.require_pass_distance and "distance_filter" in coloc:
        keep &= coloc["distance_filter"].eq("PASS")
    if args.set_label:
        keep &= coloc["set_label"].fillna("").str.contains(args.set_label, case=False, regex=False)

    targets = coloc.loc[keep].copy()
    targets = targets.dropna(subset=["QTL_Beta", "GWAS_SNP_Beta", "PP.H4.abf", "gene_name"])
    targets = targets[(targets["QTL_Beta"] != 0) & (targets["GWAS_SNP_Beta"] != 0)]

    targets["risk_expression_direction"] = np.sign(targets["QTL_Beta"] * targets["GWAS_SNP_Beta"])
    targets["protective_expression_direction"] = -targets["risk_expression_direction"]
    targets["protective_direction_label"] = targets["protective_expression_direction"].map(sign_label)
    targets["risk_direction_label"] = targets["risk_expression_direction"].map(sign_label)
    targets["target_weight"] = targets["PP.H4.abf"].clip(lower=0) * (
        targets["GWAS_SNP_Beta"].abs() / targets["GWAS_SNP_SE"].replace(0, np.nan)
    ).fillna(1.0)

    sort_cols = ["PP.H4.abf", "gene_name", "GWAS", "locus"]
    targets = targets.sort_values(sort_cols, ascending=[False, True, True, True])

    out_cols = [
        "disease",
        "GWAS",
        "locus",
        "gene_name",
        "ensembl_id",
        "set_label",
        "reference",
        "feature",
        "GWAS_SNP",
        "QTL_SNP",
        "GWAS_SNP_Beta",
        "GWAS_SNP_SE",
        "QTL_Beta",
        "QTL_SE",
        "PP.H4.abf",
        "protective_expression_direction",
        "protective_direction_label",
        "risk_expression_direction",
        "risk_direction_label",
        "target_weight",
        "distance_filter",
    ]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    targets[out_cols].to_csv(args.out, sep="\t", index=False)

    gene_summary = (
        targets.assign(weighted_direction=targets["protective_expression_direction"] * targets["target_weight"])
        .groupby("gene_name", as_index=False)
        .agg(
            n_coloc=("PP.H4.abf", "size"),
            max_h4=("PP.H4.abf", "max"),
            mean_h4=("PP.H4.abf", "mean"),
            protective_score=("weighted_direction", "sum"),
            total_weight=("target_weight", "sum"),
        )
    )
    gene_summary["protective_expression_direction"] = np.sign(gene_summary["protective_score"])
    gene_summary["protective_direction_label"] = gene_summary["protective_expression_direction"].map(sign_label)
    gene_summary["abs_protective_score"] = gene_summary["protective_score"].abs()
    gene_summary = gene_summary.sort_values(["max_h4", "abs_protective_score"], ascending=False)
    gene_summary.to_csv(args.summary_out, sep="\t", index=False)

    print(f"Wrote {len(targets):,} coloc target rows to {args.out}")
    print(f"Wrote {len(gene_summary):,} gene-level targets to {args.summary_out}")


if __name__ == "__main__":
    main()
