#!/usr/bin/env python3
"""Fit/apply landmark-to-all-gene imputation for CPA fallback outputs."""

from __future__ import annotations

import argparse

import _bootstrap  # noqa: F401
from lincs_microglial.landmark_impute import apply_landmark_model, fit_ridge_landmark_to_all


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    fit = sub.add_parser("fit")
    fit.add_argument("--out-npz", default="data/processed/cpa/landmark_to_all_ridge.npz")
    fit.add_argument("--max-signatures", type=int, default=50000)
    fit.add_argument("--alpha", type=float, default=100.0)
    fit.add_argument("--seed", type=int, default=1)
    fit.add_argument("--chunk-size", type=int, default=1024)

    apply = sub.add_parser("apply")
    apply.add_argument("--landmark-pred-h5ad", required=True)
    apply.add_argument("--model-npz", default="data/processed/cpa/landmark_to_all_ridge.npz")
    apply.add_argument("--out-h5ad", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "fit":
        fit_ridge_landmark_to_all(
            out_npz=args.out_npz,
            max_signatures=args.max_signatures,
            alpha=args.alpha,
            seed=args.seed,
            chunk_size=args.chunk_size,
        )
        print(f"Wrote landmark imputation model: {args.out_npz}")
    elif args.command == "apply":
        apply_landmark_model(args.landmark_pred_h5ad, args.model_npz, args.out_h5ad)
        print(f"Wrote all-gene imputed predictions: {args.out_h5ad}")


if __name__ == "__main__":
    main()
