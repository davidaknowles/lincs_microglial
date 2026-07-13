#!/usr/bin/env bash
#SBATCH -J cpa_val_pilot
#SBATCH -p gpu
#SBATCH --gres=gpu:b6k:1
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH -o logs/cpa_validation_pilot.%j.out
#SBATCH -e logs/cpa_validation_pilot.%j.err

set -euo pipefail
module load Python/3.10.8-GCCcore-12.2.0
. "${HOME}/venv/cpa_lincs/bin/activate"
export PYTHONPATH="${PWD}:${PWD}/cpa_compat:${PYTHONPATH:-}"
mkdir -p logs data/processed/cpa results/cpa

python scripts/prepare_cpa_lincs.py \
  --mode validation \
  --gene-mode pilot \
  --out-h5ad data/processed/cpa/cpa_lincs_validation_pilot.h5ad \
  --out-eligible data/processed/cpa/cpa_lincs_eligible_drugs.tsv

python scripts/train_cpa_lincs.py \
  --h5ad data/processed/cpa/cpa_lincs_validation_pilot.h5ad \
  --out-dir results/cpa/validation_pilot_model \
  --max-epochs 300 \
  --batch-size 512

python scripts/predict_cpa_thp1.py \
  --model-dir results/cpa/validation_pilot_model \
  --h5ad data/processed/cpa/cpa_lincs_validation_pilot.h5ad \
  --out-h5ad data/processed/cpa/cpa_thp1_validation_pilot_predictions.h5ad \
  --out-target-long data/processed/cpa/cpa_thp1_validation_pilot_target_gene_summary.tsv \
  --response-source validated_holdout_prediction

python scripts/summarize_cpa_results.py \
  --observed-gene data/processed/lincs_thp1_drug_gene_summary.tsv \
  --predicted-gene data/processed/cpa/cpa_thp1_validation_pilot_target_gene_summary.tsv \
  --out data/processed/cpa/cpa_validation_pilot_metrics.tsv
