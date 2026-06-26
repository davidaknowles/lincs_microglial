#!/usr/bin/env python3
"""Rank LINCS compounds by weighted no-intercept MR~LINCS regression."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

from annotate_prioritized_drugs import MANUAL_ANNOTATIONS, read_repurposing_drugs


def normal_cdf(z: float) -> float:
    return 0.5 * math.erfc(-z / math.sqrt(2.0))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--drug-gene", default="data/processed/lincs_thp1_protective_drug_gene_scores.tsv")
    parser.add_argument("--annotations", default="data/processed/prelim_top_lincs_thp1_protective_drugs_annotated.tsv")
    parser.add_argument("--repurposing-drugs", default="data/external/repurposing_hub/repurposing_drugs_20200324.txt")
    parser.add_argument("--out-drug", default="data/processed/lincs_thp1_mr_lincs_regression_drug_scores.tsv")
    parser.add_argument("--out-gene", default="data/processed/lincs_thp1_mr_lincs_regression_gene_data.tsv")
    return parser.parse_args()


def weighted_no_intercept(frame: pd.DataFrame) -> pd.Series:
    x = frame["mean_z"].to_numpy(float)
    y = frame["mr_ivw_beta"].to_numpy(float)
    se = frame["mr_ivw_se"].to_numpy(float)
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(se) & (se > 0) & (x != 0)
    if ok.sum() < 2:
        return pd.Series(
            {
                "n_genes": int(ok.sum()),
                "slope": np.nan,
                "slope_se": np.nan,
                "slope_z": np.nan,
                "p_protective": np.nan,
                "p_two_sided": np.nan,
                "sum_wx2": np.nan,
            }
        )
    x = x[ok]
    y = y[ok]
    se = se[ok]
    w = 1.0 / np.square(se)
    sum_wx2 = np.sum(w * np.square(x))
    slope = np.sum(w * x * y) / sum_wx2
    slope_se = math.sqrt(1.0 / sum_wx2)
    slope_z = slope / slope_se
    # Protective alternative is slope < 0: drug-induced expression moves opposite to AD-risk-increasing expression.
    p_protective = normal_cdf(slope_z)
    p_two_sided = math.erfc(abs(slope_z) / math.sqrt(2.0))
    return pd.Series(
        {
            "n_genes": int(ok.sum()),
            "slope": slope,
            "slope_se": slope_se,
            "slope_z": slope_z,
            "p_protective": p_protective,
            "p_two_sided": p_two_sided,
            "sum_wx2": sum_wx2,
        }
    )


def annotate(drug: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
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
        if mask.any():
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
    return drug


def main() -> None:
    args = parse_args()
    gene = pd.read_csv(args.drug_gene, sep="\t")
    gene = gene.dropna(subset=["mean_z", "mr_ivw_beta", "mr_ivw_se"]).copy()
    gene = gene[gene["mr_ivw_se"].gt(0)].copy()
    gene["mr_weight"] = 1.0 / np.square(gene["mr_ivw_se"])
    gene["risk_product"] = gene["mean_z"] * gene["mr_ivw_beta"]
    gene["protective_product"] = -gene["risk_product"]

    drug = (
        gene.groupby(["pert_id", "pert_iname"], group_keys=False)
        .apply(weighted_no_intercept)
        .reset_index()
        .sort_values(["p_protective", "slope_z"], ascending=[True, True])
    )
    drug["protective_slope"] = -drug["slope"]
    drug["minus_log10_p_protective"] = -np.log10(drug["p_protective"].replace(0, np.nextafter(0, 1)))
    drug = annotate(drug, args)

    Path(args.out_drug).parent.mkdir(parents=True, exist_ok=True)
    drug.to_csv(args.out_drug, sep="\t", index=False)
    gene.to_csv(args.out_gene, sep="\t", index=False)
    print(f"Wrote {len(drug):,} regression drug scores to {args.out_drug}")
    print(f"Wrote {len(gene):,} regression gene data to {args.out_gene}")


if __name__ == "__main__":
    main()
