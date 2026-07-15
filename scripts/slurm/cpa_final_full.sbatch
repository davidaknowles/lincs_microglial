#!/usr/bin/env bash
#SBATCH -J cpa_final_full
#SBATCH -p gpu
#SBATCH --gres=gpu:b6k:1
#SBATCH -c 12
#SBATCH --mem=220G
#SBATCH -t 48:00:00
#SBATCH -o logs/cpa_final_full.%j.out
#SBATCH -e logs/cpa_final_full.%j.err

set -euo pipefail
module load Python/3.12.3-GCCcore-13.3.0
. "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${PWD}:${PYTHONPATH:-}"
mkdir -p logs data/processed/cpa results/cpa

python scripts/prepare_cpa_lincs.py \
  --mode final \
  --gene-mode all \
  --out-h5ad data/processed/cpa/cpa_lincs_final_all_genes.h5ad \
  --out-eligible data/processed/cpa/cpa_lincs_final_all_genes_eligible_drugs.tsv

python scripts/train_cpa_lincs.py \
  --h5ad data/processed/cpa/cpa_lincs_final_all_genes.h5ad \
  --out-dir results/cpa/final_all_genes_model \
  --max-epochs 2000 \
  --batch-size 512

python scripts/predict_cpa_thp1.py \
  --model-dir results/cpa/final_all_genes_model \
  --h5ad data/processed/cpa/cpa_lincs_final_all_genes.h5ad \
  --out-h5ad data/processed/cpa/cpa_thp1_final_all_genes_predictions.h5ad \
  --out-target-long data/processed/cpa/cpa_thp1_predicted_target_gene_summary.tsv

python scripts/rank_cpa_drugs_by_isomiga.py
