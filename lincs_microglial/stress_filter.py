"""Reusable LINCS stress/toxicity screening utilities."""

from __future__ import annotations

import h5py
import numpy as np
import pandas as pd

from lincs_microglial.lincs_io import gctx_index


STRESS_MODULES = {
    "p53_dna_damage": [
        "BAX", "BBC3", "BTG2", "CDKN1A", "DDB2", "FAS", "GADD45A", "GADD45B",
        "MDM2", "PMAIP1", "SESN1", "SESN2", "XPC",
    ],
    "apoptosis": [
        "APAF1", "BAK1", "BAX", "BBC3", "BCL2L11", "CASP3", "CASP7", "FAS",
        "PMAIP1", "TNFRSF10B",
    ],
    "integrated_stress_upr": [
        "ASNS", "ATF3", "ATF4", "CHAC1", "DDIT3", "HERPUD1", "HSPA5", "PPP1R15A",
        "TRIB3", "XBP1",
    ],
    "heat_shock": ["BAG3", "DNAJB1", "HSPA1A", "HSPA1B", "HSPB1", "HSPH1", "HSP90AA1"],
    "oxidative_stress": ["FTH1", "GCLC", "GCLM", "HMOX1", "NQO1", "SLC7A11", "SRXN1", "TXNRD1"],
    "cell_cycle_arrest": ["BTG1", "BTG2", "CDKN1A", "CDKN1B", "GADD45A", "GADD45B", "RBL2"],
}


def condition_consistency(
    target_long: pd.DataFrame,
    targets: pd.DataFrame,
    min_signatures: int,
    min_positive_fraction: float,
) -> pd.DataFrame:
    """Assess protective-direction consistency across dose/time signatures."""
    direction = targets[["gene_name", "protective_expression_direction"]].drop_duplicates("gene_name")
    long = target_long.merge(direction, on="gene_name", how="inner")
    long["protective_push_z"] = long["zscore"] * long["protective_expression_direction"]
    long["sig_protective"] = long["protective_push_z"].gt(0) & long["zscore"].abs().ge(1)
    long["sig_opposing"] = long["protective_push_z"].lt(0) & long["zscore"].abs().ge(1)
    sig_cols = [
        column
        for column in [
            "sig_id", "pert_id", "pert_iname", "pert_dose", "pert_dose_unit",
            "pert_time", "pert_time_unit",
        ]
        if column in long.columns
    ]
    signature = long.groupby(sig_cols, dropna=False, as_index=False).agg(
        n_sig_condition_protective=("sig_protective", "sum"),
        n_sig_condition_opposing=("sig_opposing", "sum"),
    )
    signature["condition_net"] = (
        signature["n_sig_condition_protective"] - signature["n_sig_condition_opposing"]
    )
    signature["condition_net_positive"] = signature["condition_net"].gt(0)
    signature["pert_dose"] = pd.to_numeric(signature.get("pert_dose"), errors="coerce")
    dose_group = [
        column
        for column in ["pert_id", "pert_time", "pert_time_unit", "pert_dose_unit"]
        if column in signature
    ]
    signature["is_lowest_dose"] = signature["pert_dose"].eq(
        signature.groupby(dose_group, dropna=False)["pert_dose"].transform("min")
    )
    n_dose_levels = signature.groupby("pert_id")["pert_dose"].nunique(dropna=True)
    lowest_positive = signature[signature["is_lowest_dose"]].groupby("pert_id")[
        "condition_net_positive"
    ].any()
    out = signature.groupby("pert_id", as_index=False).agg(
        n_thp1_signatures=("sig_id", "nunique"),
        n_positive_conditions=("condition_net_positive", "sum"),
        frac_positive_conditions=("condition_net_positive", "mean"),
        median_condition_net=("condition_net", "median"),
    )
    out["n_dose_levels"] = out["pert_id"].map(n_dose_levels).fillna(0).astype(int)
    out["lowest_dose_positive"] = out["pert_id"].map(lowest_positive).fillna(False).astype(bool)
    replicated = out["n_thp1_signatures"].ge(min_signatures)
    consistent = out["frac_positive_conditions"].ge(min_positive_fraction)
    lower_dose = out["n_dose_levels"].le(1) | out["lowest_dose_positive"]
    out["condition_consistency_assessed"] = replicated
    out["condition_consistency_pass"] = ~replicated | (consistent & lower_dose)
    out["condition_evidence"] = np.where(
        ~replicated,
        "unreplicated_single_profile",
        np.where(consistent & lower_dose, "replicated_consistent", "replicated_inconsistent"),
    )
    return out


def read_drug_mean_and_stress(
    gctx_path: str,
    gene_info_path: str,
    signature_meta: pd.DataFrame,
    chunk_size: int,
    stress_modules: dict[str, list[str]] | None = None,
) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    """Stream selected signatures and return drug means plus stress-module scores."""
    modules = STRESS_MODULES if stress_modules is None else stress_modules
    gene_info = pd.read_csv(gene_info_path, sep="\t", dtype=str)
    gene_symbols = dict(zip(gene_info["gene_id"], gene_info["gene_symbol"]))
    with h5py.File(gctx_path, "r") as h5:
        idx = gctx_index(h5)
        if idx.gene_axis != "col":
            raise ValueError("Expected a GCTX matrix with signatures stored by rows")
        sig_lookup = {sig_id: position for position, sig_id in enumerate(idx.sig_ids)}
        signature_meta = signature_meta[signature_meta["sig_id"].isin(sig_lookup)].copy()
        signature_meta["gctx_position"] = signature_meta["sig_id"].map(sig_lookup)
        signature_meta = signature_meta.sort_values("gctx_position").reset_index(drop=True)
        positions = signature_meta["gctx_position"].to_numpy(dtype=int)
        pert_ids = sorted(signature_meta["pert_id"].unique())
        pert_lookup = {pert_id: position for position, pert_id in enumerate(pert_ids)}
        pert_codes = signature_meta["pert_id"].map(pert_lookup).to_numpy(dtype=int)
        gene_names = [gene_symbols.get(gene_id, gene_id) for gene_id in idx.gene_ids]
        gene_lookup = {gene_name: position for position, gene_name in enumerate(gene_names)}
        module_indices = {
            name: np.array([gene_lookup[gene] for gene in genes if gene in gene_lookup], dtype=int)
            for name, genes in modules.items()
        }
        if any(len(indices) < 4 for indices in module_indices.values()):
            sizes = {name: len(indices) for name, indices in module_indices.items()}
            raise ValueError(f"Insufficient genes for stress modules: {sizes}")
        sums = np.zeros((len(pert_ids), len(gene_names)), dtype=np.float32)
        counts = np.bincount(pert_codes, minlength=len(pert_ids)).astype(np.float32)
        module_values = {name: np.empty(len(signature_meta), dtype=np.float32) for name in module_indices}
        for start in range(0, len(positions), chunk_size):
            stop = min(start + chunk_size, len(positions))
            block = np.asarray(idx.matrix[positions[start:stop], :], dtype=np.float32)
            np.add.at(sums, pert_codes[start:stop], block)
            for name, indices in module_indices.items():
                module_values[name][start:stop] = block[:, indices].mean(axis=1)
        drug_mean = sums / counts[:, None]
    signature_stress = signature_meta[["sig_id", "pert_id"]].copy()
    for name, values in module_values.items():
        signature_stress[f"stress_{name}"] = values
    return signature_stress, drug_mean, pert_ids


def cosine_similarity_to_centroid(matrix: np.ndarray, centroid: np.ndarray) -> np.ndarray:
    """Calculate row-wise cosine similarity to a reference centroid."""
    numerator = matrix @ centroid
    denominator = np.linalg.norm(matrix, axis=1) * np.linalg.norm(centroid)
    return np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0)
