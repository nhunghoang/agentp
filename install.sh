#!/bin/bash
# install.sh - Complete installation script for dev_agentp
# Works on macOS ARM, macOS Intel, Linux
# Usage: bash install.sh

set -euo pipefail

echo "============================================"
echo "Starting dev_agentp Installation"
echo "============================================"

INSTALL_DIR="$(pwd)"
echo "Installing in: $INSTALL_DIR"

# ----------------------------------------------------------
# [1/8] Clone or update repository
# ----------------------------------------------------------
echo ""
echo "[1/8] Setting up repository..."
cd "$INSTALL_DIR"
REPO_URL="https://github.com/nhunghoang/dev_agentp.git"
REPO_DIR="dev_agentp"
BRANCH="main"

if [ -d "$REPO_DIR/.git" ]; then
  echo "[update] $REPO_DIR"
  cd "$REPO_DIR"
  git remote set-url origin "$REPO_URL"
  git fetch --all --tags --prune
  git checkout "$BRANCH"
  git reset --hard "origin/$BRANCH"
  git submodule sync --recursive
  git submodule update --init --recursive --depth 1
else
  echo "[clone] $REPO_DIR"
  git clone --recurse-submodules --depth 1 -b "$BRANCH" "$REPO_URL" "$REPO_DIR"
  cd "$REPO_DIR"
  git submodule update --init --recursive --depth 1
fi

echo "Current commit: $(git rev-parse --short HEAD)"

# ----------------------------------------------------------
# [2/8] Install agentp Python package
# ----------------------------------------------------------
echo ""
echo "[2/8] Installing agentp package..."
cd "$INSTALL_DIR/dev_agentp"
pip uninstall -y agentp 2>/dev/null || true
pip install -e .

python - <<'PY'
import agentp, pathlib
print("agentp will import from:", pathlib.Path(agentp.__file__).parent)
PY

# ----------------------------------------------------------
# [3/8] Download TWAS models
# ----------------------------------------------------------
echo ""
echo "[3/8] Downloading TWAS models..."

PKG_DIR="$INSTALL_DIR/dev_agentp/agentp"
mkdir -p "$PKG_DIR/models"

URL="https://vanderbilt.box.com/shared/static/10hz24rk7z9r6oh7h3st84vicv6ksqfq.gz"
OUT="$PKG_DIR/models/models.tar.gz"

if command -v curl &> /dev/null; then
  curl -L "$URL" -o "$OUT"
else
  wget -O "$OUT" "$URL"
fi

tar -xzf "$OUT" --strip-components=1 -C "$PKG_DIR/models"
rm -f "$OUT"

ls -lh "$PKG_DIR/models/JTI/weights" 2>/dev/null || echo "Model weights directory missing!"

# ----------------------------------------------------------
# [4/8] Install Python requirements
# ----------------------------------------------------------
echo ""
echo "[4/8] Installing Python packages..."
pip install -r "$INSTALL_DIR/dev_agentp/requirements.txt"

# ----------------------------------------------------------
# [5/8] Install PLINK2
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

wget -q "$PLINK_URL" -O "$PLINK_DIR/plink2.zip"
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
# [6/8] Install BGENIX (macOS ARM-safe version)
# ----------------------------------------------------------
echo ""
echo "[6/8] Installing BGENIX..."

ARCH=$(uname -m)

if [[ "$OSTYPE" == "darwin"* && "$ARCH" == "arm64" ]]; then
  echo "Detected macOS ARM — using Intel (Rosetta) environment for bgenix."

  export CONDA_SUBDIR=osx-64

  if command -v conda &> /dev/null; then
    conda create -y -n bgen -c conda-forge -c bioconda bgenix
    ln -sf "$HOME/miniforge3/envs/bgen/bin/bgenix" "$HOME/.local/bin/bgenix"
  else
    echo "Error: conda not found. Please install Miniforge:"
    echo "https://github.com/conda-forge/miniforge"
    exit 1
  fi

else
  echo "Installing bgenix for architecture: $ARCH"
  conda create -y -n bgen -c conda-forge -c bioconda bgenix
  ln -sf "$HOME/miniforge3/envs/bgen/bin/bgenix" "$HOME/.local/bin/bgenix"
fi

echo "BGENIX installed:"
"$HOME/.local/bin/bgenix" -help 2>&1 | head -n5 || echo "Warning: bgenix may need verification."

# ----------------------------------------------------------
# [7/8] Install REGENIE
# ----------------------------------------------------------
echo ""
echo "[7/8] Installing REGENIE..."

if command -v conda &> /dev/null; then
  conda create -y -n gwas -c conda-forge -c bioconda regenie
  ln -sf "$HOME/miniforge3/envs/gwas/bin/regenie" "$HOME/.local/bin/regenie"
else
  echo "Error: conda not found for REGENIE."
fi

echo "REGENIE installed:"
"$HOME/.local/bin/regenie" --help 2>&1 | head -n5 || true

# ----------------------------------------------------------
# [8/8] Install R Packages
# ----------------------------------------------------------
echo ""
echo "[8/8] Installing R packages..."

if ! command -v R &> /dev/null; then
  echo "R is not installed. Please install R manually."
else
  R -q -e "install.packages(c('optparse','RColorBrewer'), repos='https://cloud.r-project.org')"
fi

echo ""
echo "============================================"
echo "Installation Complete!"
echo "============================================"
echo "Location: $INSTALL_DIR/dev_agentp"
echo ""
echo "Next steps:"
echo "Run: bash run_example.sh"
