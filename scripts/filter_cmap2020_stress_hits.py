#!/usr/bin/env python3
"""Filter protective THP1 hits for condition consistency and stress/toxicity."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import _bootstrap  # noqa: F401
from lincs_microglial.stress_filter import (
    condition_consistency,
    cosine_similarity_to_centroid,
    read_drug_mean_and_stress,
)

CYTOTOXIC_MOA_PATTERN = "|".join(
    [
        "topoisomerase inhibitor",
        "proteasome inhibitor",
        "hsp inhibitor",
        "rna polymerase inhibitor",
        "ribonucleotide reductase inhibitor",
        "thymidylate synthase inhibitor",
        "dihydrofolate reductase inhibitor",
        "dna alkylating agent",
        "dna crosslinking agent",
        "dna inhibitor",
        "rna synthesis inhibitor",
        "protein synthesis inhibitor",
        "microtubule inhibitor",
        "tubulin inhibitor",
        "tubulin polymerization inhibitor",
        "telomerase inhibitor",
        "mdm inhibitor",
        "aurora kinase inhibitor",
        "plk inhibitor",
        "bcl inhibitor",
        "atp synthase inhibitor",
    ]
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gctx",
        default="data/raw/lincs_cmap2020/level5_beta_trt_cp_n720216x12328.gctx",
    )
    parser.add_argument("--gene-info", default="data/raw/lincs_cmap2020/geneinfo_beta.txt")
    parser.add_argument(
        "--target-long",
        default="data/processed/cmap2020_thp1_target_gene_zscores.tsv.gz",
    )
    parser.add_argument("--targets", default="data/processed/protective_expression_gene_summary.tsv")
    parser.add_argument(
        "--drug-scores",
        default="data/processed/cpa/isomiga_cmap2020_cpa_abs1_protective_count_combined_drug_scores.tsv",
    )
    parser.add_argument("--min-signatures", type=int, default=2)
    parser.add_argument("--min-positive-fraction", type=float, default=0.6)
    parser.add_argument("--reference-quantile", type=float, default=0.95)
    parser.add_argument("--chunk-size", type=int, default=128)
    parser.add_argument(
        "--out-qc",
        default="data/processed/cpa/isomiga_cmap2020_abs1_stress_qc_drug_scores.tsv",
    )
    parser.add_argument(
        "--out-positive",
        default="docs/cmap2020_stress_filtered_positive_hits.tsv",
    )
    parser.add_argument(
        "--out-filtered",
        default="data/processed/cpa/isomiga_cmap2020_abs1_stress_filtered_drug_scores.tsv",
    )
    return parser.parse_args()


def failure_reason(row: pd.Series) -> str:
    reasons = []
    if not row["condition_consistency_pass"]:
        reasons.append("inconsistent_or_high_dose_only")
    if row["annotated_cytotoxic_moa"]:
        reasons.append("annotated_cytotoxic_moa")
    if not row["stress_module_pass"]:
        reasons.append("high_stress_module_score")
    if not row["cytotoxic_similarity_pass"]:
        reasons.append("cytotoxic_signature_similarity")
    return "pass" if not reasons else "|".join(reasons)


def main() -> None:
    args = parse_args()
    target_long = pd.read_csv(args.target_long, sep="\t")
    targets = pd.read_csv(args.targets, sep="\t")
    scores = pd.read_csv(args.drug_scores, sep="\t")
    observed = scores[scores["response_source"].eq("observed_thp1")].copy()
    observed["original_rank"] = np.arange(1, len(observed) + 1)

    consistency = condition_consistency(
        target_long,
        targets,
        min_signatures=args.min_signatures,
        min_positive_fraction=args.min_positive_fraction,
    )
    signature_meta = target_long.drop_duplicates("sig_id")[["sig_id", "pert_id"]]
    signature_stress, drug_mean, pert_ids = read_drug_mean_and_stress(
        args.gctx,
        args.gene_info,
        signature_meta,
        args.chunk_size,
    )

    stress_columns = [column for column in signature_stress if column.startswith("stress_")]
    drug_stress = signature_stress.groupby("pert_id", as_index=False)[stress_columns].mean()
    drug_stress["max_stress_module_score"] = drug_stress[stress_columns].max(axis=1)
    signature_stress["max_signature_stress_score"] = signature_stress[stress_columns].max(axis=1)
    stress_fraction = (
        signature_stress.assign(stress_signature=signature_stress["max_signature_stress_score"].ge(1))
        .groupby("pert_id")["stress_signature"]
        .mean()
    )
    drug_stress["frac_stress_signatures"] = drug_stress["pert_id"].map(stress_fraction)

    cytotoxic_ids = set(
        observed.loc[
            observed["moa"].fillna("").str.contains(CYTOTOXIC_MOA_PATTERN, case=False, regex=True),
            "pert_id",
        ]
    )
    pert_position = {pert_id: position for position, pert_id in enumerate(pert_ids)}
    cytotoxic_positions = [pert_position[pert_id] for pert_id in cytotoxic_ids if pert_id in pert_position]
    if len(cytotoxic_positions) < 10:
        raise ValueError(f"Only {len(cytotoxic_positions)} annotated cytotoxic reference drugs were found")
    cytotoxic_centroid = drug_mean[cytotoxic_positions].mean(axis=0)
    similarity = cosine_similarity_to_centroid(drug_mean, cytotoxic_centroid)
    similarity_by_id = dict(zip(pert_ids, similarity))

    qc = observed.merge(consistency, on="pert_id", how="left").merge(drug_stress, on="pert_id", how="left")
    qc["annotated_cytotoxic_moa"] = qc["pert_id"].isin(cytotoxic_ids)
    qc["cytotoxic_signature_similarity"] = qc["pert_id"].map(similarity_by_id)
    launched_reference = qc["clinical_phase"].eq("Launched") & ~qc["annotated_cytotoxic_moa"]
    if launched_reference.sum() < 20:
        raise ValueError(f"Only {launched_reference.sum()} launched non-cytotoxic reference drugs were found")
    stress_cutoff = qc.loc[launched_reference, "max_stress_module_score"].quantile(args.reference_quantile)
    similarity_cutoff = qc.loc[launched_reference, "cytotoxic_signature_similarity"].quantile(
        args.reference_quantile
    )
    qc["stress_module_cutoff"] = stress_cutoff
    qc["cytotoxic_similarity_cutoff"] = similarity_cutoff
    qc["stress_module_pass"] = qc["max_stress_module_score"].le(stress_cutoff)
    qc["cytotoxic_similarity_pass"] = qc["cytotoxic_signature_similarity"].le(similarity_cutoff)
    qc["stress_toxicity_filter_pass"] = (
        qc["condition_consistency_pass"]
        & ~qc["annotated_cytotoxic_moa"]
        & qc["stress_module_pass"]
        & qc["cytotoxic_similarity_pass"]
    )
    qc["stress_toxicity_filter_reason"] = qc.apply(failure_reason, axis=1)

    out_qc = Path(args.out_qc)
    out_qc.parent.mkdir(parents=True, exist_ok=True)
    qc.to_csv(out_qc, sep="\t", index=False)
    positive = qc[qc["net_sig_protective_genes"].gt(0)].copy()
    positive = positive.sort_values("original_rank")
    out_positive = Path(args.out_positive)
    out_positive.parent.mkdir(parents=True, exist_ok=True)
    positive.to_csv(out_positive, sep="\t", index=False)

    passed = positive[positive["stress_toxicity_filter_pass"]].copy()
    out_filtered = Path(args.out_filtered)
    out_filtered.parent.mkdir(parents=True, exist_ok=True)
    passed.to_csv(out_filtered, sep="\t", index=False)
    print(f"Cytotoxic reference drugs: {len(cytotoxic_positions)}")
    print(f"Launched non-cytotoxic calibration drugs: {launched_reference.sum()}")
    print(f"Stress module cutoff ({args.reference_quantile:g} quantile): {stress_cutoff:.4f}")
    print(f"Cytotoxic similarity cutoff ({args.reference_quantile:g} quantile): {similarity_cutoff:.4f}")
    print(f"Positive-net hits before filtering: {len(positive)}")
    print(f"Positive-net hits after filtering: {len(passed)}")
    print(f"Wrote {out_qc}")
    print(f"Wrote {out_positive}")
    print(f"Wrote {out_filtered}")


if __name__ == "__main__":
    main()
