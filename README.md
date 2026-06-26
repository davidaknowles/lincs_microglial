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

## Current run summary

- High-H4 expression coloc rows: 52.
- Gene-level targets: 27.
- LINCS-covered target genes: 13.
- THP1 compound Level 5 signatures: 648 signatures across 371 compound names.
- Scored compounds with at least two represented target genes: 376 compound IDs.
- Top-ranked compound by the current metric: BIBR-1532.
- Selected biology-aware examples plotted in detail: BIBR-1532, AS-605240, myriocin, AZD-8055, and triptolide.
