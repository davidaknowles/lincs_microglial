#!/usr/bin/env bash
set -euo pipefail

module load Python/3.12.3-GCCcore-13.3.0
. .venv/bin/activate

python scripts/derive_coloc_targets.py \
  --coloc isomiga_AD_coloc.txt \
  --h4-threshold 0.8 \
  --require-pass-distance \
  --set-label expression

python scripts/download_lincs_geo.py --release GSE92742

python scripts/extract_lincs_thp1_targets.py

python scripts/score_lincs_drugs.py --min-genes 2

python scripts/annotate_prioritized_drugs.py

python scripts/score_mr_lincs_drugs.py
