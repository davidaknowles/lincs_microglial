#!/usr/bin/env bash
#SBATCH -J cpa_no_top2k
#SBATCH -p gpu
#SBATCH --gres=gpu:b6k:1
#SBATCH -c 12
#SBATCH --mem=160G
#SBATCH -t 48:00:00
#SBATCH -o logs/cpa_nosmiles_top2000.%j.out
#SBATCH -e logs/cpa_nosmiles_top2000.%j.err

set -euo pipefail
module load Python/3.12.3-GCCcore-13.3.0
. "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${PWD}:${PYTHONPATH:-}"
mkdir -p logs data/processed/cpa results/cpa

python scripts/prepare_cpa_lincs.py \
  --mode final \
  --gene-mode landmark \
  --time-hours all \
  --require-smiles false \
  --drug-selection top-n \
  --top-n-drugs 2000 \
  --targets none \
  --out-h5ad data/processed/cpa/cpa_nosmiles_top2000_final_landmark.h5ad \
  --out-eligible data/processed/cpa/cpa_nosmiles_top2000_final_eligible_drugs.tsv

python scripts/train_cpa_lincs.py \
  --h5ad data/processed/cpa/cpa_nosmiles_top2000_final_landmark.h5ad \
  --out-dir results/cpa/nosmiles_top2000_final_model \
  --max-epochs 2000 \
  --batch-size 512 \
  --use-rdkit-embeddings false

python scripts/predict_cpa_thp1.py \
  --model-dir results/cpa/nosmiles_top2000_final_model \
  --h5ad data/processed/cpa/cpa_nosmiles_top2000_final_landmark.h5ad \
  --out-h5ad data/processed/cpa/cpa_nosmiles_top2000_final_landmark_predictions.h5ad \
  --out-target-long data/processed/cpa/cpa_nosmiles_top2000_final_landmark_target_gene_summary.tsv \
  --use-rdkit-embeddings false \
  --response-source nosmiles_top2000_final_landmark_prediction

python scripts/impute_l1000_to_all_genes.py fit \
  --out-npz data/processed/cpa/landmark_to_all_ridge.npz

python scripts/impute_l1000_to_all_genes.py apply \
  --landmark-pred-h5ad data/processed/cpa/cpa_nosmiles_top2000_final_landmark_predictions.h5ad \
  --model-npz data/processed/cpa/landmark_to_all_ridge.npz \
  --out-h5ad data/processed/cpa/cpa_nosmiles_top2000_final_all_gene_predictions.h5ad

python scripts/export_cpa_predictions_long.py \
  --h5ad data/processed/cpa/cpa_nosmiles_top2000_final_all_gene_predictions.h5ad \
  --targets data/processed/protective_expression_gene_summary.tsv \
  --out data/processed/cpa/cpa_nosmiles_top2000_final_target_gene_summary.tsv \
  --response-source nosmiles_top2000_final_all_gene_prediction

python scripts/rank_cpa_drugs_by_isomiga.py \
  --predicted-gene data/processed/cpa/cpa_nosmiles_top2000_final_target_gene_summary.tsv \
  --out-prefix data/processed/cpa/isomiga_cpa_nosmiles_top2000_protective_count
