# ISOMIGA AD coloc x LINCS THP1 preliminary drug matching

This directory contains a reproducible preliminary analysis for matching ISOMIGA microglia AD-colocalizing eQTLs to LINCS L1000 THP1 compound perturbations. The goal is to identify THP1 drug perturbations that move AD-colocalizing microglial eQTL genes in the genetically protective direction.

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
5. `scripts/annotate_prioritized_drugs.py`
   - Adds known MOA/target annotations from Broad Drug Repurposing Hub where available.
   - Adds manual annotations for high-priority LINCS Phase I compounds missing from the Repurposing Hub table.
   - Writes selected-drug match plots for compounds with a strong score and/or compelling microglial biology.
6. `scripts/score_mr_lincs_drugs.py`
   - Ranks drugs using MR beta magnitude and MR standard error, not only genetic direction.
   - Reports `mr_lincs_effect = -LINCS_mean_z * mr_ivw_beta`, where positive values indicate predicted reduction in AD risk.
   - Reports `mr_lincs_precision_match = -LINCS_mean_z * (mr_ivw_beta / mr_ivw_se)`, which upweights genes with more precise naive MR estimates.
7. `scripts/score_mr_lincs_regression.py`
   - Fits a heteroskedastic no-intercept regression for each drug: `MR_beta_g = slope_drug * LINCS_z_drug,g + error_g`.
   - Uses weights `1 / MR_se_g^2`.
   - Ranks by the one-sided significance of a protective negative slope, `p(slope < 0)`.

## CMap2020 THP1 expansion

The current observed analysis uses the CMap2020 compound release rather than relying only on GSE92742. CMap2020 contains 15,814 THP1 Level 5 compound signatures across 1,667 perturbation IDs. The extraction averages all available THP1 signatures for each perturbation-gene pair, as in the original observed workflow.

Download and extract CMap2020:

```bash
python scripts/download_lincs_cmap2020.py
python scripts/extract_lincs_thp1_targets.py \
  --gctx data/raw/lincs_cmap2020/level5_beta_trt_cp_n720216x12328.gctx \
  --sig-info data/raw/lincs_cmap2020/siginfo_beta.txt \
  --gene-info data/raw/lincs_cmap2020/geneinfo_beta.txt \
  --out-long data/processed/cmap2020_thp1_target_gene_zscores.tsv.gz \
  --out-summary data/processed/cmap2020_thp1_drug_gene_summary.tsv \
  --out-coverage data/processed/cmap2020_lincs_target_gene_coverage.tsv
python scripts/combine_lincs_observed_releases.py
```

`combine_lincs_observed_releases.py` prefers CMap2020 for overlapping perturbation IDs and retains one GSE92742-only observed drug. Existing CPA predictions are then used only for IDs absent from both observed releases; CPA was not retrained for this update.

The `abs(mean_z) >= 1` combined ranking contains 2,938 perturbation IDs: 1,668 observed and 1,270 CPA-predicted. Of these, 159 observed drugs have more thresholded protective than opposing gene effects. No CPA-predicted effect reaches the threshold.

The current thresholded outputs are:

- `data/processed/cpa/isomiga_cmap2020_cpa_abs1_protective_count_combined_drug_scores.tsv`
- `docs/cmap2020_top30_biology_review.tsv`: top-30 score table with threshold-driving genes and a biological reasonableness assessment for every hit.
- `results/figures/cmap2020_cpa_top30_abs1_coloc_gene_drug_heatmap.pdf`
- `results/figures/cmap2020_cpa_top30_abs1_coloc_gene_drug_heatmap.png`

The review labels hits as `plausible_with_caveats`, `likely_nonspecific_or_toxic`, or `uninterpretable`. These labels are qualitative triage based on known pharmacology, relevance to myeloid/neurodegenerative biology, and whether the expression pattern is consistent with broad stress or cytotoxicity. They are not evidence of efficacy.

### Stress/toxicity triage

The recommended observed-hit table applies a simple post-ranking screen intended to remove signatures likely dominated by acute cell stress or cytotoxicity:

```bash
python scripts/filter_cmap2020_stress_hits.py
python scripts/plot_cpa_coloc_gene_drug_heatmap.py \
  --drug-scores data/processed/cpa/isomiga_cmap2020_abs1_stress_filtered_drug_scores.tsv \
  --gene-scores data/processed/cpa/isomiga_cmap2020_cpa_abs1_protective_count_combined_drug_gene_scores.tsv \
  --out-pdf results/figures/cmap2020_stress_filtered_top30_abs1_coloc_gene_drug_heatmap.pdf \
  --out-png results/figures/cmap2020_stress_filtered_top30_abs1_coloc_gene_drug_heatmap.png \
  --plot-data data/processed/cpa/cmap2020_stress_filtered_top30_abs1_heatmap_data.tsv \
  --min-abs-z 1
python scripts/review_cmap2020_top_hits.py \
  --drug-scores data/processed/cpa/isomiga_cmap2020_abs1_stress_filtered_drug_scores.tsv \
  --out docs/cmap2020_stress_filtered_top30_biology_review.tsv
```

The screen uses four checks:

1. Exclude annotated mechanisms whose intended acute action is cytotoxic or cytostatic, including DNA/RNA/protein synthesis disruption, proteasome or HSP inhibition, mitotic inhibition, and direct p53/apoptosis activation.
2. Calculate six stress modules from the full 12,328-gene THP1 signatures: p53/DNA damage, apoptosis, integrated stress/UPR, heat shock, oxidative stress, and cell-cycle arrest.
3. Calculate cosine similarity between each compound's full-transcriptome mean response and a centroid formed from annotated cytotoxic compounds.
4. For compounds with at least two THP1 profiles, require at least 60% of dose/time conditions to have a positive protective-minus-opposing gene count and require support at the lowest tested dose. Single-profile compounds remain eligible but are labeled `unreplicated_single_profile` rather than being treated as consistent.

Stress-module and cytotoxic-centroid thresholds are the 95th percentiles among launched compounds without an annotated cytotoxic MOA. This makes the thresholds relative to clinically used reference compounds rather than arbitrary z-score cutoffs. LINCS TAS is not used as a toxicity measure. The screen reduced 159 positive-net observed hits to 67; it is a triage filter, not a substitute for viability, morphology, or primary microglial dose-response experiments.

Outputs:

- `docs/cmap2020_stress_filtered_positive_hits.tsv`: all 159 positive-net hits, filter metrics, pass/fail status, and explicit failure reasons.
- `data/processed/cpa/isomiga_cmap2020_abs1_stress_filtered_drug_scores.tsv`: 67 passing compounds in the original genetics/LINCS rank order.
- `docs/cmap2020_stress_filtered_top30_biology_review.tsv`: biological discussion and evidence level for every filtered top-30 hit.
- `results/figures/cmap2020_stress_filtered_top30_abs1_coloc_gene_drug_heatmap.pdf`: filtered top-30 coloc-gene heatmap.

Heatmap colors use a symmetric `asinh(LINCS z)` mapping, while the legend remains labeled in the original z-score units. This preserves extreme values without compressing the more common effects around `|z| = 1-3`.

## CPA THP1 imputation

The CPA extension predicts THP1 drug responses for LINCS compounds that are well represented in other cell lines. It is kept separate from the preliminary observed-only pipeline because CPA has a different Python dependency stack.

Set up the CPA environment:

```bash
# Expects a Python-3.12-compatible CPA checkout at ../CPA by default.
# Set CPA_SOURCE=/path/to/CPA to use a different checkout.
scripts/setup_cpa_env.sh
```

Run the validation pilot:

```bash
sbatch scripts/slurm/cpa_validation_pilot.sbatch
```

Run the original SMILES/RDKit final all-gene model:

```bash
sbatch scripts/slurm/cpa_final_full.sbatch
```

If the all-gene model is not feasible after lowering batch size, run the landmark fallback:

```bash
sbatch scripts/slurm/cpa_final_landmark_fallback.sbatch
```

Run the no-SMILES learned-embedding workflow modeled after the original CPA L1000 analysis:

```bash
sbatch scripts/slurm/cpa_nosmiles_smoke.sbatch
sbatch scripts/slurm/cpa_nosmiles_top1000_reproduce.sbatch
sbatch scripts/slurm/cpa_nosmiles_top2000_final_allgenes.sbatch
```

CPA data preparation uses LINCS Level 5 z-scores, `trt_cp` compound profiles, and vehicle controls recoded as `DMSO`. The original SMILES workflow uses 6-hour profiles and requires usable canonical SMILES, at least 20 non-THP1 signatures, and at least 3 non-THP1 cell lines. The no-SMILES workflow uses learned perturbation embeddings instead of RDKit embeddings, selects drugs by non-THP1 profile count, and therefore can include compounds without SMILES. Synthetic THP1 query rows use THP1 vehicle-control expression as the basal state, with the target drug and a nonnegative raw-dose `cpa_dose` encoded for CPA; `log_dose` is retained as metadata.

The original CPA L1000 paper analysis used 978 measured landmark genes and learned drug embeddings. The local no-SMILES reproduction follows that design with the local GSE92742 release, where the available treatment metadata has 71 cell lines rather than the paper notebook's 82-cell-line prepared object.

Two split modes are implemented:

- Validation mode holds out all observed THP1 compound perturbations as `ood`, while retaining THP1 controls. This tests whether CPA can recover measured THP1 responses from non-THP1 drug evidence.
- Final mode trains with all observed THP1 perturbations and predicts only eligible drugs without observed THP1 compound profiles.

The default final attempt trains CPA on all 12,328 LINCS genes. If this is infeasible, the fallback trains CPA on measured L1000 landmark genes and then applies a ridge landmark-to-all-gene imputation model fit from observed LINCS signatures. Fallback outputs label genes as landmark-measured or post-CPA-imputed.

The no-SMILES final workflow always trains CPA on 978 measured landmark genes first, then imputes the THP1 predictions to all 12,328 LINCS genes with the same ridge landmark-to-all-gene model before ISOMIGA ranking.

After CPA prediction, `scripts/rank_cpa_drugs_by_isomiga.py` combines observed and predicted THP1 responses with ISOMIGA coloc targets. Observed THP1 responses are preferred when both observed and predicted responses exist. Drugs are ranked by the simple count of matched ISOMIGA genes moved in the protective direction:

```text
protective_push_z = mean_z * protective_expression_direction
n_protective_genes = count(protective_push_z > 0)
```

The combined ranking writes observed-only, predicted-only, and combined drug tables under `data/processed/cpa/`.

Earlier GSE92742/CPA tables retained for comparison:

- `data/processed/cpa/isomiga_cpa_nosmiles_top2000_protective_count_combined_drug_scores.tsv`: primary combined ranking. Observed THP1 responses are used when available; CPA-imputed THP1 responses fill in drugs without observed THP1 profiles.
- `data/processed/cpa/isomiga_cpa_nosmiles_top2000_protective_count_observed_only_drug_scores.tsv`: observed THP1 drugs only.
- `data/processed/cpa/isomiga_cpa_nosmiles_top2000_protective_count_predicted_only_drug_scores.tsv`: CPA-imputed THP1 drugs only.
- `docs/cpa_nosmiles_top30_biology_review.tsv`: concise biology review for the top 30 primary ranked drugs.
- `results/figures/cpa_nosmiles_top30_coloc_gene_drug_heatmap.pdf`: heatmap of top ranked drugs by coloc genes. Rows are sorted by naive MR `mr_ivw_beta`; tile color is the drug effect on expression (`mean_z`), and black dots mark gene-drug effects in the genetically protective direction.
- `data/processed/cpa/isomiga_cpa_nosmiles_top2000_abs1_protective_count_combined_drug_scores.tsv`: sensitivity ranking requiring `abs(mean_z) >= 1` for significant protective/opposing gene counts. The primary rank is `net_sig_protective_genes = n_sig_protective_genes - n_sig_opposing_genes`.
- `results/figures/cpa_nosmiles_top30_abs1_coloc_gene_drug_heatmap.pdf`: heatmap from the `abs(mean_z) >= 1` ranking. Tiles still show continuous expression z-scores; ticks mark protective effects passing the absolute z-score threshold, and crosses mark thresholded opposing effects.

The earlier GSE92742 top 30 and biology review are retained as historical outputs. They should not be treated as the current ranking now that measured CMap2020 THP1 responses are available.

The `abs(mean_z) >= 1` sensitivity ranking changes the observed-drug order and promotes drugs with more thresholded protective than opposing gene effects. Predicted-only CPA drugs do not enter the top 30 under this threshold because their all-gene imputed ISOMIGA target effects are shrunken: the maximum absolute predicted target-gene `mean_z` is below 1 in the current top-2,000 no-SMILES run.

## Gene-drug matching logic

### Genetic protective direction

For each high-colocalization eQTL row, the pipeline compares the AD GWAS SNP beta and eQTL SNP beta for the coloc SNP:

- `QTL_Beta > 0` means the allele increases expression of the target gene.
- `GWAS_SNP_Beta > 0` means the allele increases AD risk.
- If `QTL_Beta` and `GWAS_SNP_Beta` have the same sign, increased expression tracks higher AD risk, so decreased expression is inferred as protective.
- If the signs differ, increased expression tracks lower AD risk, so increased expression is inferred as protective.

In code:

```text
risk_expression_direction = sign(QTL_Beta * GWAS_SNP_Beta)
protective_expression_direction = -risk_expression_direction
```

Rows are filtered to gene-expression colocalizations with `PP.H4.abf >= 0.8` and `distance_filter == PASS`. Multiple coloc rows for the same gene are collapsed to a gene-level target by summing the protective direction across rows:

```text
target_weight = 1
gene protective_score = sum(protective_expression_direction * target_weight)
```

The sign of the gene-level `protective_score` is the target direction used for drug matching. `PP.H4.abf` is used as a high-confidence filter and reporting column, not as a drug-matching weight, because it varies only from 0.8 to 1.0 after filtering.

### Naive Mendelian randomization estimate

The pipeline also reports a very naive single-variant MR/Wald ratio for each coloc row:

```text
mr_wald_beta = GWAS_SNP_Beta / QTL_Beta
```

Using first-order delta-method error propagation:

```text
mr_wald_se = sqrt(
    (GWAS_SNP_SE / QTL_Beta)^2
    + (GWAS_SNP_Beta * QTL_SE / QTL_Beta^2)^2
)
```

At the gene level, coloc-row Wald estimates are combined with a fixed-effect inverse-variance weighted estimate:

```text
mr_ivw_beta = sum(mr_wald_beta / mr_wald_se^2) / sum(1 / mr_wald_se^2)
mr_ivw_se = sqrt(1 / sum(1 / mr_wald_se^2))
```

Interpretation is deliberately cautious. This is not a full MR analysis: it reuses coloc summary rows, may duplicate correlated evidence across references/GWAS releases, and does not model pleiotropy or LD among instruments. It is included as an effect-size sanity check: positive `mr_*_beta` means genetically increased expression is associated with increased AD risk, while negative means genetically increased expression is associated with lower AD risk.

### LINCS drug direction

LINCS Level 5 signatures are moderated z-scores. For each THP1 compound and represented target gene, the pipeline averages all THP1 `trt_cp` signatures for that compound-gene pair:

```text
mean_z = mean(Level 5 z-score across THP1 signatures for compound and gene)
```

The drug-gene protective push is:

```text
protective_push_z = mean_z * protective_expression_direction
```

So `protective_push_z > 0` means the drug moves that gene in the genetically protective direction. A value `>= 1` is counted as a strong protective push.

### Compound score

Each compound is summarized across LINCS-covered coloc target genes:

- `n_target_genes`: number of represented coloc target genes.
- `n_protective_genes`: genes with `protective_push_z > 0`.
- `n_strong_protective_genes`: genes with `protective_push_z >= 1`.
- `weighted_mean_protective_push_z`: currently equal-gene-weight mean of `protective_push_z`; the column name is retained for compatibility with earlier outputs.
- `fraction_genes_protective`: `n_protective_genes / n_target_genes`.

The default ranking emphasizes multi-gene support: first `n_strong_protective_genes`, then `weighted_mean_protective_push_z`, then `fraction_genes_protective`.

This ranking is intentionally preliminary. A drug with several strong gene pushes can still have a negative mean if it strongly opposes other target genes. Those cases are useful pathway comparators but should not be treated as clean protective candidates.

### MR-aware compound score

The MR-aware ranking uses the magnitude and uncertainty of the naive gene-level MR estimates:

```text
mr_lincs_effect = -LINCS_mean_z * mr_ivw_beta
mr_lincs_precision_match = -LINCS_mean_z * (mr_ivw_beta / mr_ivw_se)
```

Positive values mean that the LINCS drug perturbation changes expression in a direction predicted by the naive MR estimate to reduce AD risk. The default MR-aware ranking sorts by `sum_precision_match`, then `mean_mr_lincs_effect`, then `fraction_positive_mr_effect`.

This score is distinct from the direction-only protective-push score. It can promote drugs whose effects fall on genes with larger or more precise MR estimates, even if their direction-only rank is lower.

### Heteroskedastic no-intercept regression score

The regression score treats the LINCS drug signature as the predictor and the naive MR beta as the outcome across target genes:

```text
MR_beta_g = slope_drug * LINCS_z_drug,g + error_g
weight_g = 1 / MR_se_g^2
```

No intercept is fit because the null drug signature (`LINCS_z = 0`) should imply no genetically predicted risk effect. A protective drug should have a negative slope: genes increased by the drug tend to have negative MR betas, and genes decreased by the drug tend to have positive MR betas.

The weighted no-intercept estimate is:

```text
slope = sum(weight_g * LINCS_z_g * MR_beta_g) / sum(weight_g * LINCS_z_g^2)
slope_se = sqrt(1 / sum(weight_g * LINCS_z_g^2))
```

Drugs are ranked by the one-sided p-value for `slope < 0`, with smaller `p_protective` indicating a more significant protective-direction match.

## MOA/target annotations

`data/processed/prelim_top_lincs_thp1_protective_drugs_annotated.tsv` adds:

- `moa`
- `target`
- `clinical_phase`
- `disease_area`
- `indication`
- `annotation_source`
- `biology_note`

Most annotations come from Broad Drug Repurposing Hub files in `data/external/repurposing_hub`, which provide mechanism-of-action and protein target annotations. Some top LINCS Phase I compounds do not appear in that Repurposing Hub release or have incomplete annotations; for those, `scripts/annotate_prioritized_drugs.py` applies a small manual fallback for high-interest compounds such as BIBR-1532, AS-605240, triptolide, myriocin, AZD-8055, and glibenclamide.

Manual notes are intended for prioritization context, not clinical interpretation.

## Main outputs

- `data/processed/protective_expression_targets.tsv`: high-H4 coloc target rows with protective direction calls.
- `data/processed/protective_expression_gene_summary.tsv`: gene-level target summary.
- `data/processed/lincs_target_gene_coverage.tsv`: which target genes are present in LINCS.
- `data/processed/lincs_thp1_protective_drug_scores.tsv`: ranked compound scores.
- `data/processed/lincs_thp1_protective_drug_gene_scores.tsv`: per-compound, per-gene protective-push scores.
- `data/processed/prelim_top_lincs_thp1_protective_drugs.tsv`: top 50 table exported by the notebook.
- `data/processed/prelim_top_lincs_thp1_protective_drugs_annotated.tsv`: top 50 table with known MOA/targets and biology notes.
- `data/processed/lincs_thp1_mr_lincs_drug_scores.tsv`: MR-aware drug ranking.
- `data/processed/lincs_thp1_mr_lincs_drug_gene_scores.tsv`: per-drug, per-gene MR-aware match scores.
- `data/processed/lincs_thp1_mr_lincs_regression_drug_scores.tsv`: weighted no-intercept regression drug ranking.
- `data/processed/lincs_thp1_mr_lincs_regression_gene_data.tsv`: per-drug, per-gene data used by the regression ranking.
- `notebooks/isomiga_lincs_prelim_analysis.ipynb`: plotnine analysis notebook.
- `notebooks/isomiga_lincs_prelim_analysis.executed.ipynb`: executed copy.
- `results/figures/*.png`: preliminary figures.

Key figures:

- `target_direction_counts.png`: counts of high-H4 genes where increased vs decreased expression is genetically protective.
- `top_coloc_targets.png`: top coloc targets by PPF.H4.
- `top_lincs_protective_drugs.png`: ranked compound-level protective push.
- `top_drug_gene_protective_push_heatmap.png`: per-gene match for the top compounds.
- `selected_promising_drug_gene_match_bars.png`: per-gene protective push for selected interesting drugs.
- `selected_promising_drug_genetics_lincs_scatter.png`: genetic target weight versus LINCS expression effect for selected drugs.
- `selected_promising_drug_summary.png`: summary scores for selected interesting drugs.
- `top_mr_lincs_drug_scores.png`: top drugs by MR-aware LINCS match.
- `top_mr_lincs_gene_driver_scatter.png`: per-gene drivers of the top MR-aware drug matches.
- `top_mr_lincs_regression_drug_scores.png`: drugs ranked by one-sided significance of protective regression slope.
- `top_mr_lincs_regression_fits.png`: weighted no-intercept fits for the top protective slopes.

## Current run summary

- High-H4 expression coloc rows: 52.
- Gene-level targets: 27.
- LINCS-covered target genes: 13.
- CMap2020 THP1 compound Level 5 signatures: 15,814 signatures across 1,667 compound IDs.
- Combined scored compounds with at least two represented target genes: 2,938 perturbation IDs.
- Drugs with more protective than opposing effects at `abs(mean_z) >= 1`: 159.
- Top-ranked compound by the current metric: BIBR-1532.
- The updated top-30 thresholded heatmap is `results/figures/cmap2020_cpa_top30_abs1_coloc_gene_drug_heatmap.pdf`.
