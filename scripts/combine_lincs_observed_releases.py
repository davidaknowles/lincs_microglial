#!/usr/bin/env python3
"""Combine observed THP1 releases, preferring CMap2020 for overlapping drugs."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", default="data/processed/cmap2020_thp1_drug_gene_summary.tsv")
    parser.add_argument("--fallback", default="data/processed/lincs_thp1_drug_gene_summary.tsv")
    parser.add_argument("--out", default="data/processed/cmap2020_plus_gse92742_thp1_drug_gene_summary.tsv")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    primary = pd.read_csv(args.primary, sep="\t")
    fallback = pd.read_csv(args.fallback, sep="\t")
    primary["observed_release"] = "CMap2020"
    fallback = fallback.loc[~fallback["pert_id"].isin(primary["pert_id"])].copy()
    fallback["observed_release"] = "GSE92742_fallback"
    combined = pd.concat([primary, fallback], ignore_index=True, sort=False)
    combined = combined.sort_values(["pert_iname", "gene_name", "pert_id"])
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(out, sep="\t", index=False)
    print(f"Wrote {combined['pert_id'].nunique():,} observed THP1 drugs to {out}")
    release_counts = (
        combined[["pert_id", "observed_release"]]
        .drop_duplicates()["observed_release"]
        .value_counts()
    )
    print(release_counts.to_string())


if __name__ == "__main__":
    main()
