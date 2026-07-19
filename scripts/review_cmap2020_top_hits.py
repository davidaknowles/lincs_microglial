#!/usr/bin/env python3
"""Add biological reasonableness notes to the CMap2020 top-hit table."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


REVIEWS = {
    "BRD-K04623885": (
        "likely_nonspecific_or_toxic",
        "Telomerase inhibition can induce senescence and DNA-damage programs. The clean 6:0 gene pattern is strong, "
        "but there is no direct microglial AD rationale and the THP1 response may primarily reflect cellular stress.",
    ),
    "BRD-K18190982": (
        "uninterpretable",
        "The 5:1 pattern is favorable, but neither a resolved compound identity nor MOA/target annotation is available. "
        "Treat this as a transcriptional lead requiring chemical identity and target confirmation.",
    ),
    "BRD-K89014967": (
        "plausible_with_caveats",
        "MEK-ERK signaling regulates myeloid inflammatory activation, making MEK1/2 inhibition biologically plausible. "
        "The 6:3 pattern is mixed, and systemic MEK inhibition has substantial on-target toxicity.",
    ),
    "BRD-A35588707": (
        "likely_nonspecific_or_toxic",
        "Teniposide is a cytotoxic TOP2 inhibitor. Its 4:1 pattern is more likely to reflect DNA-damage and cell-cycle "
        "stress than a selective AD-protective microglial mechanism.",
    ),
    "BRD-K01868942": (
        "plausible_with_caveats",
        "HTR4 signaling is relevant to CNS and immune physiology, but an HTR4-antagonist mechanism is not an established "
        "microglial AD strategy. The 4:1 pattern is interesting but mechanistically weak without target validation in THP1.",
    ),
    "BRD-K01779529": (
        "likely_nonspecific_or_toxic",
        "Fluoropyruvate perturbs core pyruvate metabolism. The 4:1 pattern is therefore likely driven by broad "
        "bioenergetic stress rather than a tractable, selective protective mechanism.",
    ),
    "BRD-K93095519": (
        "likely_nonspecific_or_toxic",
        "SJ-172550 disrupts MDM4-p53 signaling, but it has chemical-stability/promiscuity concerns and reported neuronal "
        "toxicity. Its 3:0 pattern is more consistent with a p53/stress response than an AD-repurposing lead.",
    ),
    "BRD-A98248982": (
        "uninterpretable",
        "The 3:0 pattern is favorable, but the perturbagen has no resolved name, MOA, or target annotation. Chemical "
        "identity and orthogonal target data are required before biological interpretation.",
    ),
    "BRD-K30677119": (
        "uninterpretable",
        "PP-30 has a 3:0 pattern but lacks a reliable public MOA annotation in the current resources. It remains a "
        "chemical-signature lead rather than a biologically interpretable candidate.",
    ),
    "BRD-K75532464": (
        "plausible_with_caveats",
        "FTI-276 inhibits protein farnesylation and can alter RAS-family signaling, vesicle trafficking, and inflammatory "
        "programs. Those processes are relevant to myeloid biology, but prenylation inhibition is broad and the 3:0 result "
        "needs confirmation with better-characterized inhibitors.",
    ),
    "BRD-K71799778": (
        "plausible_with_caveats",
        "BML-259 inhibits CDK5 and CDK2. CDK5 has neurodegeneration relevance, but a THP1 signature may also reflect CDK2 "
        "or broader kinase effects; the 3:0 pattern supports use as a pathway probe, not yet as a therapeutic lead.",
    ),
    "BRD-K16701932": (
        "likely_nonspecific_or_toxic",
        "Gly-Gly-PALO inhibits ornithine transcarbamylase and perturbs arginine/urea-cycle metabolism. The 3:0 pattern may "
        "reflect metabolic stress, with little direct support for selective microglial AD benefit.",
    ),
    "BRD-K57080016": (
        "plausible_with_caveats",
        "Selumetinib is a clinically characterized MEK1/2 inhibitor. Its 3:0 pattern and inflammatory MAPK rationale are "
        "coherent, although chronic CNS exposure, microglial specificity, and class toxicity remain concerns.",
    ),
    "BRD-K37694030": (
        "plausible_with_caveats",
        "Doxepin is CNS-penetrant and modulates histamine and monoamine signaling, but sedation and anticholinergic burden "
        "limit an AD-repurposing argument. The 3:0 THP1 pattern has no clear target-level link to the coloc genes.",
    ),
    "BRD-K73395020": (
        "plausible_with_caveats",
        "ARP-101 is an MMP2 inhibitor that can induce autophagy-associated cell death. Extracellular-matrix remodeling is "
        "relevant to neuroinflammation, but the mixed 4:2 pattern may reflect autophagic or cytotoxic stress.",
    ),
    "BRD-A73680854": (
        "plausible_with_caveats",
        "PT-630 inhibits FAP and related dipeptidyl peptidases, which can affect immune signaling. Its 4:2 pattern is mixed, "
        "and there is no strong evidence that this extracellular-protease mechanism is beneficial in microglia or AD.",
    ),
    "BRD-A13122391": (
        "plausible_with_caveats",
        "Triptolide suppresses NF-kB/inflammatory transcription, which is relevant to activated myeloid cells. However, it "
        "is a broad transcriptional inhibitor with substantial toxicity, and the 4:2 pattern includes opposing effects.",
    ),
    "BRD-A09533288": (
        "plausible_with_caveats",
        "Verapamil has a coherent calcium-signaling and vascular/immune rationale and is clinically used. The 4:2 pattern "
        "is mixed, and peripheral THP1 effects do not establish beneficial CNS or microglial exposure.",
    ),
    "BRD-K59325863": (
        "likely_nonspecific_or_toxic",
        "Delanzomib is a proteasome inhibitor. Proteostasis is relevant to neurodegeneration, but proteasome blockade is "
        "cytotoxic and the 4:2 pattern likely captures a general stress response.",
    ),
    "BRD-K15108141": (
        "likely_nonspecific_or_toxic",
        "Gemcitabine is an antiproliferative nucleoside analog targeting DNA synthesis. The 4:2 pattern is likely driven by "
        "replication/metabolic stress and provides little rationale for chronic AD treatment.",
    ),
    "BRD-K83354763": (
        "uninterpretable",
        "The 3:1 pattern is positive, but the perturbagen lacks a resolved compound name, MOA, and target annotation. It "
        "requires chemical identification before prioritization.",
    ),
    "BRD-K53903639": (
        "uninterpretable",
        "CHEMBL-1222381 has a 3:1 pattern but no usable MOA or target annotation in the current resources. The result is "
        "not biologically interpretable without identity and pharmacology curation.",
    ),
    "BRD-A36275421": (
        "uninterpretable",
        "The 3:1 pattern is positive, but the perturbagen has no resolved name, MOA, or target annotation. Treat it as an "
        "uncharacterized signature hit pending chemical and target validation.",
    ),
    "BRD-K93176058": (
        "plausible_with_caveats",
        "AC-55649 is an RAR-alpha/beta agonist. Retinoid signaling can alter myeloid differentiation and inflammatory state, "
        "making the 3:1 pattern plausible, but it is preclinical and includes an opposing USP6NL effect.",
    ),
    "BRD-K12762134": (
        "plausible_with_caveats",
        "XAV-939 inhibits tankyrase and alters Wnt/beta-catenin and PAR-dependent signaling. These pathways intersect with "
        "inflammation and glial state, but the mechanism is broad and the 3:1 pattern needs microglial validation.",
    ),
    "BRD-A19500257": (
        "likely_nonspecific_or_toxic",
        "Geldanamycin inhibits HSP90 and produces a broad proteotoxic response. Despite a 3:1 pattern, toxicity and heat-shock "
        "pathway activation make it more useful as a mechanism probe than a repurposing candidate.",
    ),
    "BRD-K55395145": (
        "likely_nonspecific_or_toxic",
        "Pemetrexed is a cytotoxic antifolate that blocks nucleotide synthesis. Its 3:1 pattern is likely an antiproliferative "
        "or metabolic-stress signature rather than selective correction of AD biology.",
    ),
    "BRD-K96123349": (
        "plausible_with_caveats",
        "Brequinar inhibits DHODH and can suppress activated immune cells through pyrimidine depletion. This immunometabolic "
        "mechanism is plausible, but the 3:1 pattern may also reflect antiproliferative stress and is not microglia-specific.",
    ),
    "BRD-K27316855": (
        "plausible_with_caveats",
        "Calcitriol activates VDR and has established immunomodulatory effects in myeloid cells. The 3:1 pattern is coherent, "
        "but systemic calcium effects, CNS exposure, and the opposing PICALM response require caution.",
    ),
    "BRD-K71726959": (
        "uninterpretable",
        "The 3:1 pattern is favorable, but no resolved compound identity, MOA, or target is available. It should remain an "
        "uncharacterized transcriptional lead until the perturbagen is identified.",
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--drug-scores",
        default="data/processed/cpa/isomiga_cmap2020_cpa_abs1_protective_count_combined_drug_scores.tsv",
    )
    parser.add_argument(
        "--gene-scores",
        default="data/processed/cpa/isomiga_cmap2020_cpa_abs1_protective_count_combined_drug_gene_scores.tsv",
    )
    parser.add_argument("--top-n", type=int, default=30)
    parser.add_argument("--out", default="docs/cmap2020_top30_biology_review.tsv")
    return parser.parse_args()


def joined_genes(frame: pd.DataFrame, column: str) -> str:
    return "|".join(frame.loc[frame[column], "gene_name"].astype(str))


def main() -> None:
    args = parse_args()
    drugs = pd.read_csv(args.drug_scores, sep="\t").head(args.top_n).copy()
    genes = pd.read_csv(args.gene_scores, sep="\t")
    missing = sorted(set(drugs["pert_id"]) - set(REVIEWS))
    if missing:
        raise ValueError(f"Missing manual reviews for top-hit perturbations: {missing}")

    drivers = []
    for pert_id, frame in genes[genes["pert_id"].isin(drugs["pert_id"])].groupby("pert_id"):
        drivers.append(
            {
                "pert_id": pert_id,
                "protective_genes_abs1": joined_genes(frame, "significant_protective"),
                "opposing_genes_abs1": joined_genes(frame, "significant_opposing"),
            }
        )
    drugs.insert(0, "rank", range(1, len(drugs) + 1))
    drugs["reasonableness"] = drugs["pert_id"].map(lambda x: REVIEWS[x][0])
    drugs["biology_read"] = drugs["pert_id"].map(lambda x: REVIEWS[x][1])
    out = drugs.merge(pd.DataFrame(drivers), on="pert_id", how="left")

    ordered = [
        "rank",
        "pert_id",
        "pert_iname",
        "response_source",
        "net_sig_protective_genes",
        "n_sig_protective_genes",
        "n_sig_opposing_genes",
        "protective_genes_abs1",
        "opposing_genes_abs1",
        "clinical_phase",
        "moa",
        "target",
        "disease_area",
        "indication",
        "reasonableness",
        "biology_read",
    ]
    out = out[ordered]
    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, sep="\t", index=False)
    print(f"Wrote {len(out)} reviewed top hits to {path}")


if __name__ == "__main__":
    main()
