# ISOMIGA AD coloc x LINCS THP1 preliminary drug matching

This directory contains a reproducible preliminary analysis for matching ISOMIGA microglia AD-colocalizing eQTLs to LINCS L1000 THP1 compound perturbations.

## Environment

Use the cluster module system and the project venv:

```bash
module load Python/3.12.3-GCCcore-13.3.0
. .venv/bin/activate
```

The venv was created from the module Python and includes pandas, numpy, scipy, h5py, plotnine, Jupyter, and tqdm.

## Pipeline

Run the whole workflow with:

```bash
scripts/run_pipeline.sh
```

Steps:

1. `scripts/derive_coloc_targets.py`
   - Filters `isomiga_AD_coloc.txt` to gene-expression rows with `PP.H4.abf >= 0.8` and `distance_filter == PASS`.
   - Infers protective expression direction from the signs of `GWAS_SNP_Beta` and `QTL_Beta`.
2. `scripts/download_lincs_geo.py --release GSE92742`
   - Downloads LINCS Phase I metadata and Level 5 imputed signatures from GEO.
   - GSE92742 is used because it contains THP1 compound signatures; GSE70138 metadata has no THP1 `trt_cp` signatures.
3. `scripts/extract_lincs_thp1_targets.py`
   - Extracts only LINCS-covered target genes from the large GCTX matrix for THP1 compound signatures.
4. `scripts/score_lincs_drugs.py`
   - Scores compounds by whether drug-induced expression pushes target genes in the inferred protective direction.

## Main outputs

- `data/processed/protective_expression_targets.tsv`: high-H4 coloc target rows with protective direction calls.
- `data/processed/protective_expression_gene_summary.tsv`: gene-level target summary.
- `data/processed/lincs_target_gene_coverage.tsv`: which target genes are present in LINCS.
- `data/processed/lincs_thp1_protective_drug_scores.tsv`: ranked compound scores.
- `data/processed/lincs_thp1_protective_drug_gene_scores.tsv`: per-compound, per-gene protective-push scores.
- `data/processed/prelim_top_lincs_thp1_protective_drugs.tsv`: top 50 table exported by the notebook.
- `notebooks/isomiga_lincs_prelim_analysis.ipynb`: plotnine analysis notebook.
- `notebooks/isomiga_lincs_prelim_analysis.executed.ipynb`: executed copy.
- `results/figures/*.png`: preliminary figures.

## Current run summary

- High-H4 expression coloc rows: 52.
- Gene-level targets: 27.
- LINCS-covered target genes: 13.
- THP1 compound Level 5 signatures: 648 signatures across 371 compound names.
- Scored compounds with at least two represented target genes: 376 compound IDs.
- Top-ranked compound by the current metric: BIBR-1532.
