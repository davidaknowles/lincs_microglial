"""Protective-count ranking for observed and CPA-predicted THP1 drug responses."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


TARGET_COLS = [
    "gene_name",
    "protective_expression_direction",
    "protective_direction_label",
    "max_h4",
    "mean_h4",
    "n_coloc",
    "total_weight",
    "max_gwas_z_abs",
    "mr_ivw_beta",
    "mr_ivw_se",
    "mr_ivw_p",
]


def load_repurposing_annotations(path: str | Path | None) -> pd.DataFrame:
    if path is None or not Path(path).exists():
        return pd.DataFrame(columns=["name_key", "clinical_phase", "moa", "target", "disease_area", "indication"])
    rep = pd.read_csv(path, sep="\t", skiprows=9, dtype=str)
    rep["name_key"] = rep["pert_iname"].str.lower()
    cols = ["name_key", "clinical_phase", "moa", "target", "disease_area", "indication"]
    return rep[[c for c in cols if c in rep.columns]].drop_duplicates("name_key")


def combine_observed_predicted(
    observed_gene: pd.DataFrame,
    predicted_gene: pd.DataFrame | None,
) -> pd.DataFrame:
    obs = observed_gene.copy()
    obs["response_source"] = "observed_thp1"
    if predicted_gene is None or predicted_gene.empty:
        return obs
    pred = predicted_gene.copy()
    if "response_source" not in pred.columns:
        pred["response_source"] = "final_cpa_unknown_prediction"
    all_gene = pd.concat([obs, pred], ignore_index=True, sort=False)
    priority = {"observed_thp1": 0, "validated_holdout_prediction": 1, "final_cpa_unknown_prediction": 2}
    all_gene["_priority"] = all_gene["response_source"].map(priority).fillna(9)
    all_gene = (
        all_gene.sort_values("_priority")
        .drop_duplicates(["pert_id", "gene_name"], keep="first")
        .drop(columns="_priority")
    )
    return all_gene


def protective_count_rank(
    targets: pd.DataFrame,
    drug_gene: pd.DataFrame,
    min_genes: int = 2,
    min_abs_z: float = 0.0,
    repurposing_drugs: str | Path | None = "data/external/repurposing_hub/repurposing_drugs_20200324.txt",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    target_cols = [c for c in TARGET_COLS if c in targets.columns]
    merged = drug_gene.merge(targets[target_cols], on="gene_name", how="inner")
    merged["protective_push_z"] = merged["mean_z"] * merged["protective_expression_direction"]
    merged["pushes_protective"] = merged["protective_push_z"].gt(0)
    merged["passes_abs_z"] = merged["mean_z"].abs().ge(min_abs_z)
    merged["significant_protective"] = merged["pushes_protective"] & merged["passes_abs_z"]
    merged["significant_opposing"] = merged["protective_push_z"].lt(0) & merged["passes_abs_z"]
    merged["sum_abs_z_protective_component"] = np.where(
        merged["pushes_protective"], merged["mean_z"].abs(), 0.0
    )
    merged["sum_abs_z_sig_protective_component"] = np.where(
        merged["significant_protective"], merged["mean_z"].abs(), 0.0
    )
    group_cols = ["pert_id", "pert_iname", "response_source"]
    drug = (
        merged.groupby(group_cols, as_index=False)
        .agg(
            n_matched_isomiga_genes=("gene_name", "nunique"),
            n_protective_genes=("pushes_protective", "sum"),
            n_opposing_genes=("pushes_protective", lambda x: int((~x).sum())),
            n_sig_protective_genes=("significant_protective", "sum"),
            n_sig_opposing_genes=("significant_opposing", "sum"),
            frac_protective_genes=("pushes_protective", "mean"),
            frac_sig_protective_genes=("significant_protective", "mean"),
            sum_abs_z_protective=("sum_abs_z_protective_component", "sum"),
            sum_abs_z_sig_protective=("sum_abs_z_sig_protective_component", "sum"),
            mean_lincs_z_matched=("mean_z", "mean"),
            min_n_signatures=("n_signatures", "min"),
        )
    )
    drug["net_sig_protective_genes"] = drug["n_sig_protective_genes"] - drug["n_sig_opposing_genes"]
    drug["min_abs_z_threshold"] = min_abs_z
    drug = drug[drug["n_matched_isomiga_genes"].ge(min_genes)].copy()
    rep = load_repurposing_annotations(repurposing_drugs)
    if not rep.empty:
        drug["name_key"] = drug["pert_iname"].str.lower()
        drug = drug.merge(rep, on="name_key", how="left").drop(columns=["name_key"])
    for col in ["moa", "target", "clinical_phase", "disease_area", "indication"]:
        if col not in drug.columns:
            drug[col] = "not annotated" if col != "indication" else ""
        else:
            drug[col] = drug[col].fillna("not annotated" if col != "indication" else "")
    drug = drug.sort_values(
        [
            "net_sig_protective_genes",
            "n_sig_protective_genes",
            "frac_sig_protective_genes",
            "sum_abs_z_sig_protective",
            "n_protective_genes",
            "frac_protective_genes",
            "sum_abs_z_protective",
            "n_matched_isomiga_genes",
        ],
        ascending=False,
    )
    merged = merged.sort_values(["pert_iname", "protective_push_z"], ascending=[True, False])
    return drug, merged


def write_rankings(
    targets_path: str | Path,
    observed_gene_path: str | Path,
    predicted_gene_path: str | Path | None,
    out_prefix: str | Path,
    min_genes: int = 2,
    min_abs_z: float = 0.0,
    repurposing_drugs: str | Path | None = "data/external/repurposing_hub/repurposing_drugs_20200324.txt",
) -> None:
    targets = pd.read_csv(targets_path, sep="\t")
    observed = pd.read_csv(observed_gene_path, sep="\t")
    predicted = pd.read_csv(predicted_gene_path, sep="\t") if predicted_gene_path and Path(predicted_gene_path).exists() else None
    combined = combine_observed_predicted(observed, predicted)
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    for label, frame in [
        ("observed_only", combined[combined["response_source"].eq("observed_thp1")]),
        ("predicted_only", combined[~combined["response_source"].eq("observed_thp1")]),
        ("combined", combined),
    ]:
        if frame.empty:
            continue
        drug, gene = protective_count_rank(
            targets,
            frame,
            min_genes=min_genes,
            min_abs_z=min_abs_z,
            repurposing_drugs=repurposing_drugs,
        )
        drug.to_csv(out_prefix.with_name(f"{out_prefix.name}_{label}_drug_scores.tsv"), sep="\t", index=False)
        gene.to_csv(out_prefix.with_name(f"{out_prefix.name}_{label}_drug_gene_scores.tsv"), sep="\t", index=False)
