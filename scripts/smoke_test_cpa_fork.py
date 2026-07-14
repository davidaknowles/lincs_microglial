#!/usr/bin/env python
"""Synthetic CPA smoke test for the fork-installed package."""

import numpy as np
import pandas as pd
import torch
import anndata as ad
import cpa


def main() -> None:
    print("torch", torch.__version__, torch.version.cuda, "cuda_available", torch.cuda.is_available(), flush=True)
    print("cpa", cpa.__file__, flush=True)
    if torch.cuda.is_available():
        print("device", torch.cuda.get_device_name(0), flush=True)

    rng = np.random.default_rng(1)
    n_obs = 240
    n_vars = 32
    conditions = np.array(["DMSO"] * 80 + ["drug_a"] * 80 + ["drug_b"] * 80)
    split = np.array(["train"] * 180 + ["test"] * 30 + ["ood"] * 30)
    rng.shuffle(split)

    obs = pd.DataFrame(
        {
            "condition_ID": conditions,
            "log_dose": np.where(conditions == "DMSO", 0.0, 1.0),
            "smiles_rdkit": np.where(
                conditions == "drug_a",
                "CCO",
                np.where(conditions == "drug_b", "CCN", ""),
            ),
            "cell_type": np.where(np.arange(n_obs) % 2 == 0, "THP1", "A549"),
            "split": split,
        }
    )

    x = rng.normal(size=(n_obs, n_vars)).astype("float32")
    x[conditions == "drug_a", :4] += 0.5
    x[conditions == "drug_b", 4:8] -= 0.5
    adata = ad.AnnData(x, obs=obs, var=pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)]))

    cpa.CPA.setup_anndata(
        adata,
        perturbation_key="condition_ID",
        control_group="DMSO",
        dosage_key="log_dose",
        smiles_key="smiles_rdkit",
        is_count_data=False,
        categorical_covariate_keys=["cell_type"],
    )
    model = cpa.CPA(
        adata,
        split_key="split",
        train_split="train",
        valid_split="test",
        test_split="ood",
        use_rdkit_embeddings=True,
        n_latent=8,
        recon_loss="gauss",
        doser_type="logsigm",
    )
    print("model initialized", flush=True)
    model.train(
        max_epochs=1,
        use_gpu=True,
        batch_size=64,
        check_val_every_n_epoch=1,
        early_stopping_patience=1,
        save_path=False,
        plan_kwargs={"lr": 1e-4, "adv_lr": 1e-4, "doser_lr": 1e-4},
    )
    print("fork blackwell synthetic one epoch train ok", flush=True)


if __name__ == "__main__":
    main()
