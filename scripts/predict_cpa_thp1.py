#!/usr/bin/env python3
"""Predict synthetic THP1 query rows from a trained CPA model."""

from __future__ import annotations

import argparse

import _bootstrap  # noqa: F401
from lincs_microglial.cpa_model import predict_query_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--out-h5ad", required=True)
    parser.add_argument("--out-target-long", default="data/processed/cpa/cpa_thp1_predicted_target_gene_summary.tsv")
    parser.add_argument("--targets", default="data/processed/protective_expression_gene_summary.tsv")
    parser.add_argument("--response-source", default="final_cpa_unknown_prediction")
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--cpu", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    predict_query_rows(
        model_dir=args.model_dir,
        h5ad=args.h5ad,
        out_h5ad=args.out_h5ad,
        out_target_long=args.out_target_long,
        targets=args.targets,
        response_source=args.response_source,
        batch_size=args.batch_size,
        use_gpu=not args.cpu,
    )
    print(f"Wrote CPA predictions: {args.out_h5ad}")
    print(f"Wrote target-gene prediction table: {args.out_target_long}")


if __name__ == "__main__":
    main()
