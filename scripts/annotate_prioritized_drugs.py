#!/usr/bin/env python3
"""Annotate prioritized LINCS compounds with MOA/targets and make match plots."""

from __future__ import annotations

import argparse
import textwrap
from pathlib import Path

import pandas as pd
from plotnine import (
    aes,
    element_text,
    facet_wrap,
    geom_col,
    geom_hline,
    geom_point,
    geom_segment,
    geom_text,
    ggplot,
    labs,
    position_dodge,
    scale_color_manual,
    scale_fill_gradient2,
    scale_fill_manual,
    theme,
    theme_bw,
)


MANUAL_ANNOTATIONS = {
    "BIBR-1532": {
        "moa": "telomerase inhibitor",
        "target": "TERT/telomerase",
        "clinical_phase": "Preclinical/tool",
        "disease_area": "oncology/tool compound",
        "annotation_source": "manual: telomerase inhibitor literature/opnMe",
        "biology_note": "Strongest transcriptomic match; tool telomerase inhibitor, so prioritize as a signature lead rather than an immediately plausible microglia therapeutic.",
    },
    "AS-605240": {
        "moa": "PI3K-gamma inhibitor",
        "target": "PIK3CG",
        "clinical_phase": "Preclinical/tool",
        "disease_area": "inflammation/immunology",
        "annotation_source": "manual: PI3K-gamma inhibitor literature",
        "biology_note": "Interesting for microglia because PI3K-gamma regulates myeloid inflammatory signaling; several strong protective-push genes despite mixed total score.",
    },
    "KU-0063794": {
        "moa": "mTOR inhibitor",
        "target": "MTOR",
        "clinical_phase": "Preclinical/tool",
        "disease_area": "oncology/autophagy",
        "annotation_source": "manual: vendor/chemical biology annotation",
        "biology_note": "Multiple strong gene-level pushes but negative weighted score; keep as a pathway comparator rather than a top protective match.",
    },
    "manumycin-a": {
        "moa": "farnesyltransferase/Ras pathway inhibitor",
        "target": "FNTA|FNTB; Ras pathway",
        "clinical_phase": "Preclinical/tool",
        "disease_area": "oncology/tool compound",
        "annotation_source": "manual: chemical biology annotation",
        "biology_note": "Several strong pushes but overall weighted opposition; useful as a Ras/isoprenylation comparator.",
    },
    "SU-11652": {
        "moa": "multi-kinase inhibitor",
        "target": "VEGFR/PDGFR/KIT-family kinases",
        "clinical_phase": "Preclinical/tool",
        "disease_area": "oncology/angiogenesis",
        "annotation_source": "manual: chemical biology annotation",
        "biology_note": "Strong pushes are offset by weighted opposing genes.",
    },
    "glibenclamide": {
        "moa": "ATP-sensitive potassium channel blocker; NLRP3 inflammasome inhibitor reported",
        "target": "ABCC8|KCNJ11; NLRP3 pathway",
        "clinical_phase": "Launched",
        "disease_area": "endocrinology/inflammation",
        "annotation_source": "manual: drug label/innate-immunity literature",
        "biology_note": "Biologically interesting for myeloid inflammasome biology, but the current THP1 target-gene score is close to neutral.",
    },
    "triptolide": {
        "moa": "XPB/TFIIH inhibitor; broad transcriptional/NF-kB suppressive effects",
        "target": "ERCC3/XPB; TFIIH complex",
        "clinical_phase": "Phase 3",
        "disease_area": "inflammation/oncology tool compound",
        "annotation_source": "manual: XPB/TFIIH target literature",
        "biology_note": "Strong multi-gene protective-direction match and anti-inflammatory biology, but likely toxicity/global transcription effects make it a mechanistic lead rather than a clean therapeutic candidate.",
    },
    "SCH-79797": {
        "moa": "PAR1 antagonist; reported antibacterial/cytotoxic activity",
        "target": "F2R/PAR1; additional off-target biology reported",
        "clinical_phase": "Preclinical/tool",
        "disease_area": "thrombosis/infectious disease tool",
        "annotation_source": "manual: chemical biology annotation",
        "biology_note": "Mixed/negative weighted match; not prioritized except as a comparator.",
    },
    "myriocin": {
        "moa": "serine palmitoyltransferase inhibitor",
        "target": "SPTLC1|SPTLC2",
        "clinical_phase": "Preclinical/tool",
        "disease_area": "sphingolipid metabolism/immunology",
        "annotation_source": "manual: sphingolipid metabolism literature",
        "biology_note": "Promising biology because sphingolipid metabolism and myeloid inflammatory state are relevant to AD; positive weighted match with broad protective direction.",
    },
    "AZD-8055": {
        "moa": "mTOR kinase inhibitor",
        "target": "MTOR",
        "clinical_phase": "Phase 1",
        "disease_area": "oncology/autophagy",
        "annotation_source": "manual: mTOR inhibitor literature/Repurposing Hub style annotation",
        "biology_note": "Positive weighted match; plausible pathway relevance through autophagy and immune metabolism, though not among the very top by strong-push count.",
    },
    "PIK-90": {
        "moa": "PI3K inhibitor",
        "target": "PIK3CA|PIK3CG|PIK3CD",
        "clinical_phase": "Preclinical/tool",
        "disease_area": "inflammation/oncology tool",
        "annotation_source": "manual: chemical biology annotation",
        "biology_note": "Relevant PI3K-family comparator with modest positive weighted match.",
    },
}


PROMISING = ["BIBR-1532", "AS-605240", "myriocin", "AZD-8055", "triptolide"]


def drug_label(name: str, moa: str, width: int = 48) -> str:
    label = f"{name} ({moa})"
    return "\n".join(textwrap.wrap(label, width=width, break_long_words=False))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", default="data/processed/lincs_thp1_protective_drug_scores.tsv")
    parser.add_argument("--gene-scores", default="data/processed/lincs_thp1_protective_drug_gene_scores.tsv")
    parser.add_argument("--targets", default="data/processed/protective_expression_gene_summary.tsv")
    parser.add_argument("--repurposing-drugs", default="data/external/repurposing_hub/repurposing_drugs_20200324.txt")
    parser.add_argument("--out-table", default="data/processed/prelim_top_lincs_thp1_protective_drugs_annotated.tsv")
    parser.add_argument("--fig-dir", default="results/figures")
    parser.add_argument("--top-n", type=int, default=50)
    return parser.parse_args()


def read_repurposing_drugs(path: str) -> pd.DataFrame:
    drugs = pd.read_csv(path, sep="\t", skiprows=9, dtype=str)
    drugs["name_key"] = drugs["pert_iname"].str.lower()
    return drugs


def annotate(scores: pd.DataFrame, repurposing: pd.DataFrame) -> pd.DataFrame:
    out = scores.copy()
    out["name_key"] = out["pert_iname"].str.lower()
    rep_cols = ["name_key", "clinical_phase", "moa", "target", "disease_area", "indication"]
    out = out.merge(repurposing[rep_cols].drop_duplicates("name_key"), on="name_key", how="left")
    out["annotation_source"] = out["moa"].notna().map({True: "Broad Repurposing Hub", False: ""})
    out["biology_note"] = ""

    for name, vals in MANUAL_ANNOTATIONS.items():
        mask = out["pert_iname"].str.lower().eq(name.lower())
        if not mask.any():
            continue
        for key, value in vals.items():
            out.loc[mask, key] = value

    out["clinical_phase"] = out["clinical_phase"].fillna("not annotated")
    out["moa"] = out["moa"].fillna("not annotated")
    out["target"] = out["target"].fillna("not annotated")
    out["disease_area"] = out["disease_area"].fillna("not annotated")
    out["indication"] = out["indication"].fillna("")
    out["annotation_source"] = out["annotation_source"].replace("", "not annotated")
    out["biology_note"] = out["biology_note"].fillna("")
    out = out.drop(columns=["name_key"])
    return out


def make_match_plots(annotated: pd.DataFrame, gene_scores: pd.DataFrame, targets: pd.DataFrame, fig_dir: Path) -> None:
    selected = [x for x in PROMISING if x in set(annotated["pert_iname"])]
    selected_meta = annotated[annotated["pert_iname"].isin(selected)].copy()
    selected_meta["pert_iname"] = pd.Categorical(selected_meta["pert_iname"], categories=selected, ordered=True)

    gene_scores = gene_scores.merge(
        selected_meta[["pert_id", "pert_iname", "moa", "target", "biology_note"]],
        on=["pert_id", "pert_iname"],
        how="inner",
    )
    gene_scores["drug_label"] = [drug_label(name, moa) for name, moa in zip(gene_scores["pert_iname"], gene_scores["moa"])]
    label_order = [drug_label(row.pert_iname, row.moa) for row in selected_meta.itertuples()]
    gene_scores["drug_label"] = pd.Categorical(gene_scores["drug_label"], categories=label_order, ordered=True)
    target_cols = [
        "gene_name",
        "protective_direction_label",
        "max_h4",
        "max_gwas_z_abs",
        "mr_ivw_beta",
        "abs_protective_score",
    ]
    gene_scores = gene_scores.merge(targets[target_cols], on="gene_name", how="left", suffixes=("", "_target"))
    gene_scores["gene_label"] = gene_scores["gene_name"].astype(str)
    gene_order = (
        targets[targets["gene_name"].isin(gene_scores["gene_name"])]
        .sort_values(["max_h4", "abs_protective_score"], ascending=False)["gene_name"]
        .tolist()
    )
    gene_scores["gene_label"] = pd.Categorical(gene_scores["gene_label"], categories=gene_order[::-1], ordered=True)

    p = (
        ggplot(gene_scores, aes("gene_label", "protective_push_z", fill="protective_push_z"))
        + geom_col()
        + geom_hline(yintercept=0, color="#333333", size=0.4)
        + facet_wrap("~drug_label", ncol=1)
        + scale_fill_gradient2(low="#9c2f2f", mid="#f7f7f7", high="#2f7d4f", midpoint=0)
        + labs(
            x="ISOMIGA AD coloc target gene",
            y="LINCS z-score in protective direction",
            fill="Protective push",
            title="",
        )
        + theme_bw()
        + theme(figure_size=(6.2, 7.0), axis_text_x=element_text(rotation=45, ha="right"), strip_text=element_text(size=8))
    )
    p.save(fig_dir / "selected_promising_drug_gene_match_bars.png", dpi=300)
    p.save(fig_dir / "selected_promising_drug_gene_match_bars.pdf")

    signed = gene_scores.assign(
        genetics_direction=gene_scores["protective_direction_label"],
        lincs_direction=gene_scores["mean_z"].map(lambda x: "increase" if x > 0 else "decrease"),
        match=gene_scores["protective_push_z"].gt(0).map({True: "matches", False: "opposes"}),
    )
    p = (
        ggplot(signed, aes("max_gwas_z_abs", "mean_z", color="match"))
        + geom_hline(yintercept=0, color="#555555", size=0.4)
        + geom_point(aes(size="max_h4"), alpha=0.9)
        + geom_text(aes(label="gene_name"), nudge_y=0.15, size=7, va="bottom", show_legend=False)
        + facet_wrap("~pert_iname", ncol=2)
        + scale_color_manual(values={"matches": "#2f7d4f", "opposes": "#9c2f2f"})
        + labs(
            x="Max coloc GWAS |z| for gene",
            y="LINCS THP1 mean z-score",
            color="Drug direction",
            size="Max PPF.H4",
            title="Genetic evidence versus drug-induced expression",
        )
        + theme_bw()
        + theme(figure_size=(10, 8))
    )
    p.save(fig_dir / "selected_promising_drug_genetics_lincs_scatter.png", dpi=300)

    summary = selected_meta[
        [
            "pert_iname",
            "n_protective_genes",
            "n_strong_protective_genes",
            "weighted_mean_protective_push_z",
            "fraction_genes_protective",
        ]
    ].copy()
    summary["pert_iname"] = pd.Categorical(summary["pert_iname"], categories=selected[::-1], ordered=True)
    p = (
        ggplot(summary, aes("weighted_mean_protective_push_z", "pert_iname"))
        + geom_segment(aes(x=0, xend="weighted_mean_protective_push_z", y="pert_iname", yend="pert_iname"), color="#777777")
        + geom_point(aes(size="n_strong_protective_genes", color="fraction_genes_protective"))
        + labs(
            x="Weighted mean protective push z",
            y="Selected compound",
            size="Strong genes",
            color="Fraction protective",
            title="Selected biologically interesting LINCS matches",
        )
        + theme_bw()
        + theme(figure_size=(7, 3.5))
    )
    p.save(fig_dir / "selected_promising_drug_summary.png", dpi=300)


def main() -> None:
    args = parse_args()
    scores = pd.read_csv(args.scores, sep="\t")
    gene_scores = pd.read_csv(args.gene_scores, sep="\t")
    targets = pd.read_csv(args.targets, sep="\t")
    repurposing = read_repurposing_drugs(args.repurposing_drugs)

    annotated = annotate(scores, repurposing).head(args.top_n).copy()
    columns = [
        "pert_id",
        "pert_iname",
        "moa",
        "target",
        "clinical_phase",
        "disease_area",
        "indication",
        "annotation_source",
        "biology_note",
        "n_target_genes",
        "n_protective_genes",
        "n_strong_protective_genes",
        "n_opposing_genes",
        "weighted_mean_protective_push_z",
        "fraction_genes_protective",
        "mean_protective_push_z",
        "median_protective_push_z",
        "mean_target_mr_ivw_beta",
        "median_target_mr_ivw_beta",
        "min_target_mr_ivw_p",
        "min_n_signatures",
        "total_gene_signatures",
    ]
    Path(args.out_table).parent.mkdir(parents=True, exist_ok=True)
    annotated[columns].to_csv(args.out_table, sep="\t", index=False)

    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)
    make_match_plots(annotated, gene_scores, targets, fig_dir)
    print(f"Wrote annotated prioritized table to {args.out_table}")
    print(f"Wrote selected-drug plots to {fig_dir}")


if __name__ == "__main__":
    main()
