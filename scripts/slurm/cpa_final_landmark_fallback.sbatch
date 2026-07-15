#!/usr/bin/env bash
#SBATCH -J cpa_final_lm
#SBATCH -p gpu
#SBATCH --gres=gpu:b6k:1
#SBATCH -c 12
#SBATCH --mem=160G
#SBATCH -t 48:00:00
#SBATCH -o logs/cpa_final_landmark.%j.out
#SBATCH -e logs/cpa_final_landmark.%j.err

set -euo pipefail
module load Python/3.12.3-GCCcore-13.3.0
. "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${PWD}:${PYTHONPATH:-}"
mkdir -p logs data/processed/cpa results/cpa

python scripts/prepare_cpa_lincs.py \
  --mode final \
  --gene-mode landmark \
  --out-h5ad data/processed/cpa/cpa_lincs_final_landmark.h5ad \
  --out-eligible data/processed/cpa/cpa_lincs_final_landmark_eligible_drugs.tsv

python scripts/train_cpa_lincs.py \
  --h5ad data/processed/cpa/cpa_lincs_final_landmark.h5ad \
  --out-dir results/cpa/final_landmark_model \
  --max-epochs 2000 \
  --batch-size 512

python scripts/predict_cpa_thp1.py \
  --model-dir results/cpa/final_landmark_model \
  --h5ad data/processed/cpa/cpa_lincs_final_landmark.h5ad \
  --out-h5ad data/processed/cpa/cpa_thp1_final_landmark_predictions.h5ad \
  --out-target-long data/processed/cpa/cpa_thp1_final_landmark_target_gene_summary.tsv

python scripts/impute_l1000_to_all_genes.py fit \
  --out-npz data/processed/cpa/landmark_to_all_ridge.npz

python scripts/impute_l1000_to_all_genes.py apply \
  --landmark-pred-h5ad data/processed/cpa/cpa_thp1_final_landmark_predictions.h5ad \
  --model-npz data/processed/cpa/landmark_to_all_ridge.npz \
  --out-h5ad data/processed/cpa/cpa_thp1_final_landmark_to_all_predictions.h5ad

# Export target-gene long table from the all-gene fallback prediction.
python - <<'PY'
import scanpy as sc
from lincs_microglial.cpa_model import predictions_to_target_long

pred = sc.read_h5ad("data/processed/cpa/cpa_thp1_final_landmark_to_all_predictions.h5ad")
long = predictions_to_target_long(pred, targets="data/processed/protective_expression_gene_summary.tsv")
long["response_source"] = "final_cpa_unknown_prediction"
long.to_csv("data/processed/cpa/cpa_thp1_predicted_target_gene_summary.tsv", sep="\t", index=False)
PY

python scripts/rank_cpa_drugs_by_isomiga.py
