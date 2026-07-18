#!/usr/bin/env python3
"""Plot coloc genes by top CPA-ranked drugs as a LINCS expression heatmap."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from plotnine import (
    aes,
    element_text,
    geom_text,
    geom_tile,
    ggplot,
    labs,
    scale_fill_gradient2,
    theme,
    theme_bw,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--drug-scores",
        default="data/processed/cpa/isomiga_cpa_nosmiles_top2000_protective_count_combined_drug_scores.tsv",
    )
    parser.add_argument(
        "--gene-scores",
        default="data/processed/cpa/isomiga_cpa_nosmiles_top2000_protective_count_combined_drug_gene_scores.tsv",
    )
    parser.add_argument("--targets", default="data/processed/protective_expression_gene_summary.tsv")
    parser.add_argument("--out-pdf", default="results/figures/cpa_nosmiles_top30_coloc_gene_drug_heatmap.pdf")
    parser.add_argument("--out-png", default="results/figures/cpa_nosmiles_top30_coloc_gene_drug_heatmap.png")
    parser.add_argument("--plot-data", default="data/processed/cpa/cpa_nosmiles_top30_coloc_gene_drug_heatmap_data.tsv")
    parser.add_argument("--top-n", type=int, default=30)
    parser.add_argument("--min-abs-z", type=float, default=0.0)
    return parser.parse_args()


def drug_label(row: pd.Series) -> str:
    return str(row["pert_iname"])


def main() -> None:
    args = parse_args()
    drug_scores = pd.read_csv(args.drug_scores, sep="\t")
    gene_scores = pd.read_csv(args.gene_scores, sep="\t")
    targets = pd.read_csv(args.targets, sep="\t")

    top_drugs = drug_scores.head(args.top_n).copy()
    top_drugs["drug_label"] = top_drugs.apply(drug_label, axis=1)
    drug_order = top_drugs["drug_label"].tolist()

    target_cols = [
        "gene_name",
        "mr_ivw_beta",
        "mr_ivw_se",
        "mr_ivw_p",
        "protective_expression_direction",
        "protective_direction_label",
        "max_h4",
        "max_gwas_z_abs",
        "n_coloc",
    ]
    target_cols = [c for c in target_cols if c in targets.columns]

    plot_df = gene_scores.merge(
        top_drugs[["pert_id", "pert_iname", "response_source", "drug_label"]],
        on=["pert_id", "pert_iname", "response_source"],
        how="inner",
    ).merge(targets[target_cols].drop_duplicates("gene_name"), on="gene_name", how="left", suffixes=("", "_target"))

    plot_df["mr_ivw_beta"] = pd.to_numeric(plot_df["mr_ivw_beta"], errors="coerce")
    plot_df["mean_z"] = pd.to_numeric(plot_df["mean_z"], errors="coerce")
    plot_df["protective_expression_direction"] = pd.to_numeric(
        plot_df["protective_expression_direction"], errors="coerce"
    )
    plot_df["protective_push_z"] = plot_df["mean_z"] * plot_df["protective_expression_direction"]
    if args.min_abs_z > 0:
        marker_df = plot_df[plot_df["protective_push_z"].abs().ge(args.min_abs_z)].copy()
        marker_df["marker"] = np.where(marker_df["protective_push_z"].gt(0), "✓", "×")
        title_suffix = f" (|z| >= {args.min_abs_z:g})"
    else:
        marker_df = plot_df[plot_df["protective_push_z"].gt(0)].copy()
        marker_df["marker"] = "•"
        title_suffix = ""

    gene_meta = (
        plot_df[["gene_name", "mr_ivw_beta", "protective_direction_label"]]
        .drop_duplicates("gene_name")
        .sort_values("mr_ivw_beta", ascending=True)
    )
    gene_meta["gene_label"] = [f"{row.gene_name} ({row.mr_ivw_beta:+.2f})" for row in gene_meta.itertuples()]
    gene_order = gene_meta["gene_label"].tolist()
    label_map = dict(zip(gene_meta["gene_name"], gene_meta["gene_label"]))

    plot_df["gene_label"] = plot_df["gene_name"].map(label_map)
    plot_df["gene_label"] = pd.Categorical(plot_df["gene_label"], categories=gene_order, ordered=True)
    plot_df["drug_label"] = pd.Categorical(plot_df["drug_label"], categories=drug_order, ordered=True)
    marker_df["gene_label"] = marker_df["gene_name"].map(label_map)
    marker_df["gene_label"] = pd.Categorical(marker_df["gene_label"], categories=gene_order, ordered=True)
    marker_df["drug_label"] = pd.Categorical(marker_df["drug_label"], categories=drug_order, ordered=True)

    out_data = Path(args.plot_data)
    out_data.parent.mkdir(parents=True, exist_ok=True)
    plot_df.sort_values(["drug_label", "mr_ivw_beta"]).to_csv(out_data, sep="\t", index=False)

    max_abs = float(np.nanmax(np.abs(plot_df["mean_z"]))) if not plot_df.empty else 1.0
    max_abs = max(max_abs, 1.0)
    plot = (
        ggplot(plot_df, aes("drug_label", "gene_label", fill="mean_z"))
        + geom_tile(color="#f5f5f5", size=0.25)
        + geom_text(marker_df, aes(label="marker"), size=6, color="#111111")
        + scale_fill_gradient2(low="#2166ac", mid="#f7f7f7", high="#b2182b", midpoint=0, limits=(-max_abs, max_abs))
        + labs(
            x="Top ranked drug",
            y="ISOMIGA coloc gene, sorted by naive MR beta",
            fill="Drug effect\non expression\n(LINCS z)",
            title=f"Top CPA/LINCS drug effects on ISOMIGA AD coloc genes{title_suffix}",
        )
        + theme_bw()
        + theme(
            figure_size=(7.2, 3.7),
            axis_title=element_text(size=8),
            axis_text_x=element_text(rotation=45, ha="right", size=5),
            axis_text_y=element_text(size=6.5),
            legend_title=element_text(size=7),
            legend_text=element_text(size=6.5),
            plot_title=element_text(size=9),
        )
    )

    out_pdf = Path(args.out_pdf)
    out_png = Path(args.out_png)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    plot.save(out_pdf, verbose=False)
    plot.save(out_png, dpi=300, verbose=False)
    print(f"Wrote {out_pdf}")
    print(f"Wrote {out_png}")
    print(f"Wrote {out_data}")


if __name__ == "__main__":
    main()
