#!/usr/bin/env python3
"""Rank LINCS compounds by MR-aware drug/gene expression matching."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from annotate_prioritized_drugs import MANUAL_ANNOTATIONS, read_repurposing_drugs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--drug-gene", default="data/processed/lincs_thp1_protective_drug_gene_scores.tsv")
    parser.add_argument("--annotations", default="data/processed/prelim_top_lincs_thp1_protective_drugs_annotated.tsv")
    parser.add_argument("--repurposing-drugs", default="data/external/repurposing_hub/repurposing_drugs_20200324.txt")
    parser.add_argument("--out-drug", default="data/processed/lincs_thp1_mr_lincs_drug_scores.tsv")
    parser.add_argument("--out-gene", default="data/processed/lincs_thp1_mr_lincs_drug_gene_scores.tsv")
    parser.add_argument("--strong-threshold", type=float, default=2.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    gene = pd.read_csv(args.drug_gene, sep="\t")
    gene = gene.dropna(subset=["mean_z", "mr_ivw_beta", "mr_ivw_se"]).copy()
    gene = gene[gene["mr_ivw_se"].gt(0)].copy()

    # Positive values mean the drug-induced expression change predicts lower AD risk.
    gene["mr_lincs_effect"] = -gene["mean_z"] * gene["mr_ivw_beta"]
    gene["mr_lincs_effect_se_from_mr"] = gene["mean_z"].abs() * gene["mr_ivw_se"]
    gene["mr_lincs_precision_match"] = -gene["mean_z"] * (gene["mr_ivw_beta"] / gene["mr_ivw_se"])
    gene["mr_lincs_positive"] = gene["mr_lincs_effect"].gt(0)
    gene["mr_lincs_strong_precision"] = gene["mr_lincs_precision_match"].ge(args.strong_threshold)

    drug = (
        gene.groupby(["pert_id", "pert_iname"], as_index=False)
        .agg(
            n_genes=("gene_name", "nunique"),
            n_positive_mr_effect=("mr_lincs_positive", "sum"),
            n_strong_precision_match=("mr_lincs_strong_precision", "sum"),
            mean_mr_lincs_effect=("mr_lincs_effect", "mean"),
            median_mr_lincs_effect=("mr_lincs_effect", "median"),
            sum_mr_lincs_effect=("mr_lincs_effect", "sum"),
            mean_precision_match=("mr_lincs_precision_match", "mean"),
            median_precision_match=("mr_lincs_precision_match", "median"),
            sum_precision_match=("mr_lincs_precision_match", "sum"),
            mean_lincs_z=("mean_z", "mean"),
            mean_mr_ivw_beta=("mr_ivw_beta", "mean"),
            median_mr_ivw_beta=("mr_ivw_beta", "median"),
        )
    )
    drug["fraction_positive_mr_effect"] = drug["n_positive_mr_effect"] / drug["n_genes"].replace(0, np.nan)
    drug = drug.sort_values(
        ["sum_precision_match", "mean_mr_lincs_effect", "fraction_positive_mr_effect"],
        ascending=False,
    )

    rep_path = Path(args.repurposing_drugs)
    if rep_path.exists():
        rep = read_repurposing_drugs(str(rep_path))
        rep_cols = ["name_key", "clinical_phase", "moa", "target", "disease_area", "indication"]
        drug["name_key"] = drug["pert_iname"].str.lower()
        drug = drug.merge(rep[rep_cols].drop_duplicates("name_key"), on="name_key", how="left")
        drug["annotation_source"] = drug["moa"].notna().map({True: "Broad Repurposing Hub", False: ""})
        drug["biology_note"] = ""
        drug = drug.drop(columns=["name_key"])

    ann_path = Path(args.annotations)
    if ann_path.exists():
        ann = pd.read_csv(ann_path, sep="\t")
        ann_cols = ["pert_id", "moa", "target", "clinical_phase", "disease_area", "indication", "annotation_source", "biology_note"]
        ann = ann[[c for c in ann_cols if c in ann.columns]].drop_duplicates("pert_id")
        drug = drug.merge(ann, on="pert_id", how="left", suffixes=("", "_existing"))
        for col in ["moa", "target", "clinical_phase", "disease_area", "indication", "annotation_source", "biology_note"]:
            existing = f"{col}_existing"
            if existing in drug.columns:
                drug[col] = drug[col].fillna(drug[existing])
                drug = drug.drop(columns=[existing])

    for name, vals in MANUAL_ANNOTATIONS.items():
        mask = drug["pert_iname"].str.lower().eq(name.lower())
        if not mask.any():
            continue
        for key, value in vals.items():
            drug.loc[mask, key] = value

    for col in ["moa", "target", "clinical_phase", "disease_area", "annotation_source"]:
        if col not in drug:
            drug[col] = "not annotated"
        else:
            drug[col] = drug[col].fillna("not annotated").replace("", "not annotated")
    if "biology_note" not in drug:
        drug["biology_note"] = ""
    else:
        drug["biology_note"] = drug["biology_note"].fillna("")

    Path(args.out_drug).parent.mkdir(parents=True, exist_ok=True)
    drug.to_csv(args.out_drug, sep="\t", index=False)
    gene.sort_values(["pert_iname", "mr_lincs_precision_match"], ascending=[True, False]).to_csv(
        args.out_gene, sep="\t", index=False
    )
    print(f"Wrote {len(drug):,} MR-LINCS drug scores to {args.out_drug}")
    print(f"Wrote {len(gene):,} MR-LINCS drug-gene scores to {args.out_gene}")


if __name__ == "__main__":
    main()
