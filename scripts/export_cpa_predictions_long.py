#!/usr/bin/env python3
"""Export CPA prediction AnnData to long-form gene response table."""

from __future__ import annotations

import argparse
from pathlib import Path

import scanpy as sc

import _bootstrap  # noqa: F401
from lincs_microglial.cpa_model import predictions_to_target_long


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--targets", default=None, help="Optional target-gene table. Omit to export all genes.")
    parser.add_argument("--response-source", default="final_cpa_unknown_prediction")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    pred = sc.read_h5ad(args.h5ad)
    long = predictions_to_target_long(pred, targets=args.targets, response_source=args.response_source)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    long.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote long-form CPA predictions: {args.out}")


if __name__ == "__main__":
    main()
