#!/usr/bin/env bash
#SBATCH -J cpa_nosmoke
#SBATCH -p gpu
#SBATCH --gres=gpu:b6k:1
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 01:00:00
#SBATCH -o logs/cpa_nosmiles_smoke.%j.out
#SBATCH -e logs/cpa_nosmiles_smoke.%j.err

set -euo pipefail
module load Python/3.12.3-GCCcore-13.3.0
. "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${PWD}:${PYTHONPATH:-}"
mkdir -p logs data/processed/cpa results/cpa

python scripts/prepare_cpa_lincs.py \
  --mode validation \
  --gene-mode landmark \
  --time-hours all \
  --require-smiles false \
  --drug-selection top-n \
  --top-n-drugs 100 \
  --targets none \
  --split-strategy paper \
  --out-h5ad data/processed/cpa/cpa_nosmiles_smoke_top100_landmark.h5ad \
  --out-eligible data/processed/cpa/cpa_nosmiles_smoke_top100_eligible_drugs.tsv

python scripts/train_cpa_lincs.py \
  --h5ad data/processed/cpa/cpa_nosmiles_smoke_top100_landmark.h5ad \
  --out-dir results/cpa/nosmiles_smoke_top100_model \
  --max-epochs 2 \
  --batch-size 512 \
  --use-rdkit-embeddings false

python scripts/predict_cpa_thp1.py \
  --model-dir results/cpa/nosmiles_smoke_top100_model \
  --h5ad data/processed/cpa/cpa_nosmiles_smoke_top100_landmark.h5ad \
  --out-h5ad data/processed/cpa/cpa_nosmiles_smoke_top100_predictions.h5ad \
  --out-target-long data/processed/cpa/cpa_nosmiles_smoke_top100_target_gene_summary.tsv \
  --use-rdkit-embeddings false \
  --response-source nosmiles_smoke_prediction
