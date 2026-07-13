#!/usr/bin/env python3
"""Rank observed and CPA-predicted THP1 drugs by ISOMIGA protective-gene count."""

from __future__ import annotations

import argparse

import _bootstrap  # noqa: F401
from lincs_microglial.ranking import write_rankings


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", default="data/processed/protective_expression_gene_summary.tsv")
    parser.add_argument("--observed-gene", default="data/processed/lincs_thp1_drug_gene_summary.tsv")
    parser.add_argument("--predicted-gene", default="data/processed/cpa/cpa_thp1_predicted_target_gene_summary.tsv")
    parser.add_argument("--out-prefix", default="data/processed/cpa/isomiga_cpa_protective_count")
    parser.add_argument("--repurposing-drugs", default="data/external/repurposing_hub/repurposing_drugs_20200324.txt")
    parser.add_argument("--min-genes", type=int, default=2)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    write_rankings(
        targets_path=args.targets,
        observed_gene_path=args.observed_gene,
        predicted_gene_path=args.predicted_gene,
        out_prefix=args.out_prefix,
        min_genes=args.min_genes,
        repurposing_drugs=args.repurposing_drugs,
    )
    print(f"Wrote observed/predicted/combined rankings with prefix: {args.out_prefix}")


if __name__ == "__main__":
    main()
