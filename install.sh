#!/bin/bash
# install.sh - Public repo install for agentp
# Creates one conda env for python deps (agentp_env) + installs external tools.
# Usage: bash install.sh

set -euo pipefail

AGENTP_ENV="agentp_env"   # <-- hardcode python env name here

# ---- sanity ----
if ! command -v conda &> /dev/null; then
  echo "ERROR: conda not found in PATH."
  exit 1
fi

CONDA_BASE="$(conda info --base)"
mkdir -p "$HOME/.local/bin"

echo "============================================"
echo "Starting agentp Installation"
echo "============================================"

REPO_DIR="$(pwd)"
echo "Repo: $REPO_DIR"

if [[ ! -f "$REPO_DIR/pyproject.toml" && ! -f "$REPO_DIR/setup.py" ]]; then
  echo "ERROR: This doesn't look like the agentp repo (no pyproject.toml or setup.py found)."
  exit 1
fi

git -C "$REPO_DIR" rev-parse --short HEAD 2>/dev/null | awk '{print "Current commit:", $0}' || true

# ----------------------------------------------------------
# [1/8] Create / use agentp python env
# ----------------------------------------------------------
echo ""
echo "[1/8] Creating/using conda env: $AGENTP_ENV"

if conda env list | awk '{print $1}' | grep -qx "$AGENTP_ENV"; then
  echo "Env $AGENTP_ENV already exists. Skipping creation."
else
  # Choose a stable python (3.11 recommended; 3.13 has more breakage)
  conda create -y -n "$AGENTP_ENV" -c conda-forge python=3.11 pip
fi

ENV_PY="$CONDA_BASE/envs/$AGENTP_ENV/bin/python"
if [[ ! -x "$ENV_PY" ]]; then
  echo "ERROR: env python not found at: $ENV_PY"
  exit 1
fi

echo "Env python: $ENV_PY"

# ----------------------------------------------------------
# [2/8] Install Python requirements into agentp_env
# ----------------------------------------------------------
echo ""
echo "[2/8] Installing Python requirements into $AGENTP_ENV..."
conda run -n "$AGENTP_ENV" python -m pip install --upgrade pip setuptools wheel

if [[ -f "$REPO_DIR/requirements.txt" ]]; then
  conda run -n "$AGENTP_ENV" python -m pip install -r "$REPO_DIR/requirements.txt"
else
  echo "WARNING: requirements.txt not found; skipping."
fi

# ----------------------------------------------------------
# [3/8] Install agentp (editable) into agentp_env
# ----------------------------------------------------------
echo ""
echo "[3/8] Installing agentp package (editable) into $AGENTP_ENV..."
conda run -n "$AGENTP_ENV" python -m pip uninstall -y agentp 2>/dev/null || true
conda run -n "$AGENTP_ENV" python -m pip install -e "$REPO_DIR"

conda run -n "$AGENTP_ENV" python - <<PY
import os, agentp
print("agentp imported from:", os.path.realpath(agentp.__file__))
PY

# ----------------------------------------------------------
# [4/8] Download TWAS models into repo
# ----------------------------------------------------------
echo ""
echo "[4/8] Downloading TWAS models..."

PKG_DIR="$REPO_DIR/agentp"
mkdir -p "$PKG_DIR/models"

URL="https://vanderbilt.box.com/shared/static/10hz24rk7z9r6oh7h3st84vicv6ksqfq.gz"
OUT="$PKG_DIR/models/models.tar.gz"

if command -v curl &> /dev/null; then
  curl -fL "$URL" -o "$OUT"
elif command -v wget &> /dev/null; then
  wget -O "$OUT" "$URL"
else
  echo "ERROR: curl or wget required to download models."
  exit 1
fi

tar -xzf "$OUT" --strip-components=1 -C "$PKG_DIR/models"
rm -f "$OUT"

ls -lh "$PKG_DIR/models/JTI/weights" 2>/dev/null || echo "WARNING: Model weights directory missing!"

# ----------------------------------------------------------
# [5/8] Install PLINK2 into ~/.local/bin
# ----------------------------------------------------------
echo ""
echo "[5/8] Installing PLINK2..."

PLINK_DIR="$HOME/.local/bin"
mkdir -p "$PLINK_DIR"

if [[ "$OSTYPE" == "darwin"* ]]; then
  PLINK_URL="https://s3.amazonaws.com/plink2-assets/alpha6/plink2_mac_20241203.zip"
else
  PLINK_URL="https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20241203.zip"
fi

if command -v curl &> /dev/null; then
  curl -fL "$PLINK_URL" -o "$PLINK_DIR/plink2.zip"
else
  wget -q "$PLINK_URL" -O "$PLINK_DIR/plink2.zip"
fi

unzip -o -q "$PLINK_DIR/plink2.zip" -d "$PLINK_DIR"
rm "$PLINK_DIR/plink2.zip"
chmod +x "$PLINK_DIR/plink2"

# Safe PATH append (avoid permission denied)
{
  echo "export PATH=\"$PLINK_DIR:\$PATH\""
} | tee -a "$HOME/.zshrc" "$HOME/.bashrc" >/dev/null 2>&1 || true

export PATH="$PLINK_DIR:$PATH"
echo "PLINK2 installed: $($PLINK_DIR/plink2 --version | head -n1)"

# ----------------------------------------------------------
# [6/8] Install BGENIX (separate env bgen)
# ----------------------------------------------------------
echo ""
echo "[6/8] Installing BGENIX..."

ARCH=$(uname -m)
if [[ "$OSTYPE" == "darwin"* && "$ARCH" == "arm64" ]]; then
  echo "Detected macOS ARM — using Intel (Rosetta) subdir for bgenix env."
  export CONDA_SUBDIR=osx-64
fi

if ! conda env list | awk '{print $1}' | grep -qx "bgen"; then
  conda create -y -n bgen -c conda-forge -c bioconda bgenix
else
  echo "Environment bgen already exists. Skipping creation."
fi

ln -sf "$CONDA_BASE/envs/bgen/bin/bgenix" "$HOME/.local/bin/bgenix"
echo "BGENIX installed:"
"$HOME/.local/bin/bgenix" -help 2>&1 | head -n5 || echo "WARNING: bgenix may need verification."

# ----------------------------------------------------------
# [7/8] Install REGENIE (separate env gwas)
# ----------------------------------------------------------
echo ""
echo "[7/8] Installing REGENIE..."

if conda env list | awk '{print $1}' | grep -qx "gwas"; then
  echo "Environment gwas already exists. Skipping creation."
else
  if [[ -d "$CONDA_BASE/envs/gwas" && ! -f "$CONDA_BASE/envs/gwas/conda-meta/history" ]]; then
    echo "Found stale gwas directory; removing: $CONDA_BASE/envs/gwas"
    rm -rf "$CONDA_BASE/envs/gwas"
  fi
  conda create -y -n gwas -c conda-forge -c bioconda regenie
fi

ln -sf "$CONDA_BASE/envs/gwas/bin/regenie" "$HOME/.local/bin/regenie"
echo "REGENIE installed:"
"$HOME/.local/bin/regenie" --help 2>&1 | head -n5 || true

# ----------------------------------------------------------
# [8/8] Install R packages (optional)
# ----------------------------------------------------------
echo ""
echo "[8/8] Installing R packages..."
if ! command -v R &> /dev/null; then
  echo "R is not installed. Skipping."
else
  R -q -e "install.packages(c('optparse','RColorBrewer'), repos='https://cloud.r-project.org')"
fi

echo ""
echo "============================================"
echo "Installation Complete!"
echo "============================================"
echo "Repo: $REPO_DIR"
echo "Python env: $AGENTP_ENV"
echo ""
echo "Next step:"
echo "  bash run_example.sh"
