#!/usr/bin/env python3
"""Train a CPA model on a prepared LINCS AnnData file."""

from __future__ import annotations

import argparse

import _bootstrap  # noqa: F401
from lincs_microglial.cpa_model import train_cpa


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--max-epochs", type=int, default=300)
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--check-val-every-n-epoch", type=int, default=5)
    parser.add_argument("--early-stopping-patience", type=int, default=10)
    parser.add_argument("--cpu", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    train_cpa(
        h5ad=args.h5ad,
        out_dir=args.out_dir,
        max_epochs=args.max_epochs,
        batch_size=args.batch_size,
        use_gpu=not args.cpu,
        check_val_every_n_epoch=args.check_val_every_n_epoch,
        early_stopping_patience=args.early_stopping_patience,
    )
    print(f"Wrote CPA model: {args.out_dir}")


if __name__ == "__main__":
    main()
