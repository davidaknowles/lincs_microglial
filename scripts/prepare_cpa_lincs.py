#!/usr/bin/env python3
"""Build CPA-ready AnnData from LINCS GSE92742 Level 5 signatures."""

from __future__ import annotations

import argparse

import _bootstrap  # noqa: F401
from lincs_microglial.cpa_prep import prepare_cpa_data, write_anndata


def parse_bool(value: str) -> bool:
    if value.lower() in {"1", "true", "yes", "y"}:
        return True
    if value.lower() in {"0", "false", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(f"Expected boolean value, got {value!r}")


def parse_time(value: str) -> float | None:
    if value.lower() == "all":
        return None
    return float(value)


def parse_optional_path(value: str) -> str | None:
    if value.lower() in {"none", "null", ""}:
        return None
    return value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["validation", "final"], default="validation")
    parser.add_argument("--gene-mode", choices=["pilot", "all", "landmark"], default="pilot")
    parser.add_argument("--time-hours", type=parse_time, default=6)
    parser.add_argument("--min-non-thp1-signatures", type=int, default=20)
    parser.add_argument("--min-non-thp1-cell-lines", type=int, default=3)
    parser.add_argument("--require-smiles", type=parse_bool, default=True)
    parser.add_argument("--drug-selection", choices=["coverage", "top-n"], default="coverage")
    parser.add_argument("--top-n-drugs", type=int, default=None)
    parser.add_argument("--split-strategy", choices=["thp1", "paper"], default="thp1")
    parser.add_argument("--targets", type=parse_optional_path, default="data/processed/protective_expression_gene_summary.tsv")
    parser.add_argument("--out-h5ad", default="data/processed/cpa/cpa_lincs_validation_pilot.h5ad")
    parser.add_argument("--out-eligible", default="data/processed/cpa/cpa_lincs_eligible_drugs.tsv")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--chunk-size", type=int, default=2048)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    prepared = prepare_cpa_data(
        mode=args.mode,
        gene_mode=args.gene_mode,
        time_hours=args.time_hours,
        min_non_thp1_signatures=args.min_non_thp1_signatures,
        min_non_thp1_cell_lines=args.min_non_thp1_cell_lines,
        require_smiles=args.require_smiles,
        drug_selection=args.drug_selection,
        top_n_drugs=args.top_n_drugs,
        targets=args.targets,
        seed=args.seed,
        chunk_size=args.chunk_size,
        split_strategy=args.split_strategy,
    )
    write_anndata(prepared, args.out_h5ad, args.out_eligible)
    print(f"Wrote CPA AnnData: {args.out_h5ad}")
    print(f"Wrote eligible drug table: {args.out_eligible}")
    print(prepared.obs["split"].value_counts().to_string())


if __name__ == "__main__":
    main()
