#!/usr/bin/env bash
set -euo pipefail

module load Python/3.12.3-GCCcore-13.3.0

python -m venv "${HOME}/venv/cpa_blackwell"
. "${HOME}/venv/cpa_blackwell/bin/activate"
python -m pip install --upgrade pip wheel setuptools
python -m pip install torch torchaudio --index-url https://download.pytorch.org/whl/cu128
CPA_SOURCE="${CPA_SOURCE:-../CPA}"
python -m pip install -e "${CPA_SOURCE}" --extra-index-url https://download.pytorch.org/whl/cu128
python -m pip install ipykernel
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
