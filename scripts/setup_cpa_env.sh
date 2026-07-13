#!/usr/bin/env bash
set -euo pipefail

module load Python/3.10.8-GCCcore-12.2.0

python -m venv "${HOME}/venv/cpa_lincs"
. "${HOME}/venv/cpa_lincs/bin/activate"
export PYTHONPATH="${PWD}/cpa_compat:${PYTHONPATH:-}"
python -m pip install --upgrade pip wheel setuptools
python -m pip install "cpa-tools==0.8.8" ipykernel
python -m pip install "pyarrow<15"
python - <<'PY'
import anndata
import cpa
import numpy
import pandas
import scanpy
import torch

print("CPA environment ready")
print("cpa", getattr(cpa, "__version__", "unknown"))
print("torch", torch.__version__)
PY
