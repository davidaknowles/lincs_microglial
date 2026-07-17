#!/usr/bin/env bash
#SBATCH -J cpa_no_top1k
#SBATCH -p gpu
#SBATCH --gres=gpu:b6k:1
#SBATCH -c 8
#SBATCH --mem=120G
#SBATCH -t 24:00:00
#SBATCH -o logs/cpa_nosmiles_top1000.%j.out
#SBATCH -e logs/cpa_nosmiles_top1000.%j.err

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
  --top-n-drugs 1000 \
  --targets none \
  --split-strategy paper \
  --out-h5ad data/processed/cpa/cpa_nosmiles_top1000_landmark.h5ad \
  --out-eligible data/processed/cpa/cpa_nosmiles_top1000_eligible_drugs.tsv

python scripts/train_cpa_lincs.py \
  --h5ad data/processed/cpa/cpa_nosmiles_top1000_landmark.h5ad \
  --out-dir results/cpa/nosmiles_top1000_model \
  --max-epochs 500 \
  --batch-size 512 \
  --use-rdkit-embeddings false

python scripts/predict_cpa_thp1.py \
  --model-dir results/cpa/nosmiles_top1000_model \
  --h5ad data/processed/cpa/cpa_nosmiles_top1000_landmark.h5ad \
  --out-h5ad data/processed/cpa/cpa_nosmiles_top1000_predictions.h5ad \
  --out-target-long data/processed/cpa/cpa_nosmiles_top1000_target_gene_summary.tsv \
  --use-rdkit-embeddings false \
  --response-source nosmiles_top1000_prediction

python scripts/summarize_cpa_results.py \
  --observed-gene data/processed/lincs_thp1_drug_gene_summary.tsv \
  --predicted-gene data/processed/cpa/cpa_nosmiles_top1000_target_gene_summary.tsv \
  --out data/processed/cpa/cpa_nosmiles_top1000_validation_metrics.tsv

if [[ -f data/processed/cpa/cpa_thp1_validation_pilot_target_gene_summary.tsv ]]; then
  python scripts/compare_cpa_strategies.py \
    --observed-gene data/processed/lincs_thp1_drug_gene_summary.tsv \
    --a-gene data/processed/cpa/cpa_thp1_validation_pilot_target_gene_summary.tsv \
    --a-name smiles \
    --b-gene data/processed/cpa/cpa_nosmiles_top1000_target_gene_summary.tsv \
    --b-name nosmiles_top1000 \
    --out-metrics data/processed/cpa/cpa_smiles_vs_nosmiles_top1000_common_metrics.tsv \
    --out-common-gene data/processed/cpa/cpa_smiles_vs_nosmiles_top1000_common_gene.tsv
fi
