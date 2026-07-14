# Lab Notebook

## ISOMIGA x LINCS preliminary matching

The initial analysis derives AD-protective expression directions from high-H4 ISOMIGA microglial eQTL colocalizations, extracts LINCS THP1 Level 5 target-gene signatures, and ranks observed THP1 compound perturbations by the number and strength of genes moved in the protective direction.

Key components:

- `scripts/derive_coloc_targets.py`: filters coloc rows and derives gene-level protective expression directions.
- `scripts/extract_lincs_thp1_targets.py`: extracts THP1 target-gene Level 5 z-scores.
- `scripts/score_lincs_drugs.py`: direction/count-based observed THP1 drug ranking.
- `scripts/annotate_prioritized_drugs.py`: MOA/target annotation and selected plot generation.

Current finding summary:

- The observed THP1 analysis covers a subset of ISOMIGA target genes represented in LINCS.
- BIBR-1532, AS-605240, myriocin, AZD-8055, and triptolide were used as selected biology-aware examples in the preliminary plots.
- The MR and regression scripts are retained as secondary analyses, but the current CPA extension will rank drugs by simple protective-gene count.

## CPA THP1 response imputation

The CPA extension is designed to impute THP1 perturbation responses for drugs that have enough LINCS evidence in other cell lines.

Model setup:

- Input expression is LINCS GSE92742 Level 5 z-scores.
- Controls are recoded to CPA control group `DMSO`.
- The default time point is 6 hours to match THP1 conditions.
- Drug eligibility requires usable SMILES, at least 20 non-THP1 signatures, and at least 3 non-THP1 cell lines.
- CPA uses `condition_ID`, `log_dose`, `cell_type`, and RDKit embeddings from canonical SMILES.

Split strategy:

- Validation model: all observed THP1 compound perturbations are held out as `ood`; THP1 controls are retained.
- Final model: observed THP1 compound perturbations are included in training; synthetic THP1 query rows are generated for eligible drugs without observed THP1 profiles.

Gene strategy:

- First attempt: all 12,328 LINCS genes.
- If the all-gene run is infeasible, train CPA on landmark genes, then fit a ridge landmark-to-all-gene imputation model from observed LINCS signatures and apply it to CPA landmark predictions.

Ranking after CPA:

- Observed THP1 responses are preferred when available.
- CPA predictions are used for unknown THP1 drug responses.
- Drugs are ranked by `n_protective_genes`, the count of matched ISOMIGA genes with LINCS/CPA z-score sign matching the genetic protective direction.

Implementation note:

- The CPA workflow now targets a local CPA branch updated for Python 3.12, modern scvi-tools/lightning APIs, and torch CUDA 12.8 wheels for Blackwell GPUs.
- The validation pilot uses `--gene-mode pilot`, which includes L1000 landmark genes, ISOMIGA target genes, and filler genes up to 2,000 total genes. Final full mode uses all LINCS genes; the fallback mode trains on L1000 landmarks and imputes back to all genes.
- CPA training now defaults to lower optimizer learning rates (`lr`, `adv_lr`, and `doser_lr` set to `1e-4`) after the original training configuration produced unstable decoder outputs during the pilot run.
