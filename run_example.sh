# #!/bin/bash
# # run_example.sh - Example runner for agentp
# # Works on macOS ARM, macOS Intel, Linux
# # Usage: bash run_example.sh
# #
# # Assumes install.sh created a conda env named: agentp_env
# # and installed dependencies + (optionally) `pip install -e .` for this repo.

# set -euo pipefail

# echo "============================================"
# echo "Running agentp Example"
# echo "============================================"

# SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# REPO="$SCRIPT_DIR"

# # If your example lives at "$REPO/example", this is correct.
# # If your example folder name differs, change this one line.
# EX_DIR="$REPO/example"

# # Download target (only used if you want to fetch example_data)
# EX_TAR="$REPO/example.tar.gz"
# BOX_URL="https://vanderbilt.box.com/shared/static/o5s5bwwey0clt8ktabqg1abmubg2t7in.gz"

# # ------------------------------------------------------------
# # 1. AUTO-DETECT BGENIX
# # ------------------------------------------------------------
# echo ""
# echo "[Detect] Searching for bgenix..."

# BGENIX_CANDIDATES=(
#   "$HOME/miniconda3/envs/bgen/bin/bgenix"
#   "$HOME/miniforge3/envs/bgen/bin/bgenix"
#   "$HOME/miniconda3/envs/bgenix/bin/bgenix"
#   "$HOME/miniforge3/envs/bgenix/bin/bgenix"
#   "$HOME/micromamba/envs/bgen/bin/bgenix"
#   "$HOME/.local/bin/bgenix"
# )

# BGENIX_BIN=""
# for p in "${BGENIX_CANDIDATES[@]}"; do
#   if [[ -x "$p" ]]; then
#     BGENIX_BIN="$p"
#     break
#   fi
# done

# if [[ -z "$BGENIX_BIN" ]] && command -v bgenix &>/dev/null; then
#   BGENIX_BIN="$(command -v bgenix)"
# fi

# if [[ -z "$BGENIX_BIN" ]]; then
#   echo "ERROR: bgenix not found. Please run install.sh."
#   exit 1
# fi

# export PATH="$(dirname "$BGENIX_BIN"):$PATH"
# echo "[OK] Using bgenix: $BGENIX_BIN"

# # ------------------------------------------------------------
# # 2. AUTO-DETECT PLINK2
# # ------------------------------------------------------------
# echo ""
# echo "[Detect] Searching for plink2..."

# if command -v plink2 &>/dev/null; then
#   PLINK2_BIN="$(command -v plink2)"
# elif [[ -x "$HOME/.local/bin/plink2" ]]; then
#   PLINK2_BIN="$HOME/.local/bin/plink2"
#   export PATH="$HOME/.local/bin:$PATH"
# else
#   echo "ERROR: plink2 not found. Please run install.sh."
#   exit 1
# fi

# echo "[OK] Using plink2: $PLINK2_BIN"

# # ------------------------------------------------------------
# # 3. (OPTIONAL) DOWNLOAD & EXTRACT EXAMPLE DATA
# #    If example/example_data already exists, we skip download.
# # ------------------------------------------------------------
# echo ""
# if [[ -d "$EX_DIR/example_data" ]]; then
#   echo "[Data] Found existing: $EX_DIR/example_data (skip download)"
# else
#   echo "[1/4] Downloading example data..."

#   if command -v curl &>/dev/null; then
#     curl -fL "$BOX_URL" -o "$EX_TAR"
#   elif command -v wget &>/dev/null; then
#     wget -O "$EX_TAR" "$BOX_URL"
#   else
#     echo "ERROR: curl or wget required."
#     exit 1
#   fi

#   echo "Verifying archive..."
#   tar -tzf "$EX_TAR" >/dev/null

#   echo "Extracting..."
#   # Extract into repo root; archive may contain an example/ folder
#   tar -xzf "$EX_TAR" -C "$REPO"
#   rm -f "$EX_TAR"
# fi

# # ------------------------------------------------------------
# # 4. Ensure example_data path exists
# # ------------------------------------------------------------
# echo ""
# echo "[2/4] Setting up example data..."

# if [[ ! -d "$EX_DIR/example_data" ]]; then
#   echo "ERROR: example_data not found at: $EX_DIR/example_data"
#   echo "Expected your repo layout like:"
#   echo "  example/"
#   echo "    example_data/"
#   echo "    example_project/"
#   echo "    run_example.py"
#   exit 1
# fi

# echo "[OK] example_data present: $EX_DIR/example_data"

# # ------------------------------------------------------------
# # 5. CONVERT BGEN → PGEN (if needed by your pipeline)
# # ------------------------------------------------------------
# echo ""
# echo "[3/4] Converting genotype files (if needed)..."

# TMP="$EX_DIR/example_project/input_data/genotypes/tmp"
# mkdir -p "$TMP"

# DATA_DIR="$EX_DIR/example_data/genotypes"
# shopt -s nullglob

# for BGEN in "$DATA_DIR"/test_c*.bgen; do
#   base="$(basename "$BGEN" .bgen)"
#   chr="$(echo "$base" | sed -E 's/.*c([0-9]+).*/\1/')"
#   SAMPLE="$DATA_DIR/${base}.sample"
#   OUT="$TMP/c${chr}_0"

#   if [[ -f "${OUT}.pgen" && ( -f "${OUT}.pvar" || -f "${OUT}.pvar.zst" ) && -f "${OUT}.psam" ]]; then
#     echo "[skip] Already have $OUT"
#   else
#     echo "[plink2] Converting $base → $OUT"
#     plink2 --bgen "$BGEN" ref-first \
#            --sample "$SAMPLE" \
#            --make-pgen \
#            --out "$OUT"
#   fi
# done

# # ------------------------------------------------------------
# # 6. HANDLE zstd DECOMPRESSION (if needed)
# # ------------------------------------------------------------
# need_zstd=false
# for z in "$TMP"/*.pvar.zst; do
#   p="${z%.zst}"
#   if [[ -f "$z" && ! -f "$p" ]]; then need_zstd=true; fi
# done

# if $need_zstd; then
#   echo ""
#   echo "[zstd] Decompressing .pvar.zst files..."

#   if ! command -v zstd &>/dev/null && ! command -v unzstd &>/dev/null; then
#     if [[ "$OSTYPE" == "darwin"* ]]; then
#       brew install zstd
#     else
#       sudo apt-get update -qq
#       sudo apt-get install -y zstd
#     fi
#   fi

#   for z in "$TMP"/*.pvar.zst; do
#     p="${z%.zst}"
#     [[ -f "$z" && ! -f "$p" ]] && zstd -d -f "$z" -o "$p"
#   done
# fi

# echo "--- tmp contents ---"
# ls -lh "$TMP" | head || true

# # ------------------------------------------------------------
# # 7. FORCE PYTHON FROM THE SAME CONDA ENV AS install.sh
# # ------------------------------------------------------------
# echo ""
# echo "[Detect] Using agentp conda env..."

# if ! command -v conda &>/dev/null; then
#   echo "ERROR: conda not found. Please run install.sh first."
#   exit 1
# fi

# CONDA_BASE="$(conda info --base)"
# ENV_NAME="agentp_env"
# AGENTP_PY="$CONDA_BASE/envs/$ENV_NAME/bin/python"

# if [[ ! -x "$AGENTP_PY" ]]; then
#   echo "ERROR: Expected env python not found: $AGENTP_PY"
#   echo "Did you run install.sh (which creates env $ENV_NAME)?"
#   exit 1
# fi

# # Put this repo first on import path so it uses your local code
# export PYTHONPATH="$REPO${PYTHONPATH:+:$PYTHONPATH}"

# "$AGENTP_PY" - <<PY
# import sys, os
# print("python:", sys.executable)
# import agentp
# print("agentp:", os.path.realpath(agentp.__file__))
# PY

# echo "[OK] Using Python: $AGENTP_PY"

# # ------------------------------------------------------------
# # 8. RUN THE PYTHON PIPELINE
# # ------------------------------------------------------------
# echo ""
# echo "[4/4] Running example pipeline..."
# cd "$EX_DIR"

# "$AGENTP_PY" run_example.py || {
#   echo "WARNING: Pipeline returned an error. Check logs above."
#   exit 1
# }

# echo ""
# echo "============================================"
# echo "Example Run Complete!"
# echo "============================================"
# echo "Output in:"
# echo "$EX_DIR/example_project/"

#!/bin/bash
# run_example.sh - run agentp example using the env created by install.sh
# Usage: bash run_example.sh

#!/bin/bash
# run_example.sh - Run agentp example using the conda env created by install.sh
# Usage: bash run_example.sh
#
# Expects:
#   - install.sh created a conda env named "agentp_env" (editable below)
#   - install.sh installed agentp + requirements into that env
#
# This script:
#   - ensures we use the repo’s agentp (not some globally installed one)
#   - ensures example/example_data exists (download only if missing)
#   - runs example/run_example.py using the env python

set -euo pipefail

echo "============================================"
echo "Running agentp Example"
echo "============================================"

AGENTP_ENV="agentp_env"   # <-- hardcode env name (match install.sh)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$SCRIPT_DIR"

# -------------------------------
# 0) Locate conda + env python
# -------------------------------
if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found in PATH."
  echo "Please install Miniconda/Miniforge, then rerun."
  exit 1
fi

CONDA_BASE="$(conda info --base)"
ENV_PY="$CONDA_BASE/envs/$AGENTP_ENV/bin/python"

echo ""
echo "[Detect] Using agentp conda env: $AGENTP_ENV"
if [[ ! -x "$ENV_PY" ]]; then
  echo "ERROR: Expected env python not found: $ENV_PY"
  echo "Did you run install.sh (which should create env $AGENTP_ENV)?"
  exit 1
fi

# Make sure user-local tools are on PATH (plink2/bgenix/regenie symlinks live here)
export PATH="$HOME/.local/bin:$PATH"

# Force this repo to the front of import path (helps avoid importing another agentp)
export PYTHONPATH="$REPO${PYTHONPATH:+:$PYTHONPATH}"

# Quick sanity: confirm agentp is imported from this repo
echo ""
echo "[Check] Verifying agentp import location..."
conda run -n "$AGENTP_ENV" python - <<PY
import os, agentp
p = os.path.realpath(agentp.__file__)
print("agentp.__file__ =", p)
PY

# -------------------------------
# 1) Ensure example_data exists
# -------------------------------
echo ""
echo "[1/4] Setting up example data..."

EXAMPLE_DIR="$REPO/example"
EXAMPLE_DATA_DIR="$EXAMPLE_DIR/example_data"
EX_TAR="$REPO/example.tar.gz"
BOX_URL="https://vanderbilt.box.com/shared/static/o5s5bwwey0clt8ktabqg1abmubg2t7in.gz"

mkdir -p "$EXAMPLE_DIR"

if [[ -d "$EXAMPLE_DATA_DIR" ]]; then
  echo "[OK] Found existing: $EXAMPLE_DATA_DIR"
else
  echo "[Download] example_data missing, downloading archive..."
  if command -v curl >/dev/null 2>&1; then
    curl -fL "$BOX_URL" -o "$EX_TAR"
  elif command -v wget >/dev/null 2>&1; then
    wget -O "$EX_TAR" "$BOX_URL"
  else
    echo "ERROR: curl or wget required."
    exit 1
  fi

  echo "Verifying archive..."
  tar -tzf "$EX_TAR" >/dev/null

  echo "Extracting into: $REPO"
  tar -xzf "$EX_TAR" -C "$REPO"
  rm -f "$EX_TAR"

  if [[ ! -d "$EXAMPLE_DATA_DIR" ]]; then
    echo "ERROR: After extraction, still missing: $EXAMPLE_DATA_DIR"
    echo "Please check what the tarball extracted."
    exit 1
  fi
fi

# -------------------------------
# 2) Convert BGEN -> PGEN (if needed)
# -------------------------------
echo ""
echo "[2/4] Converting genotype files (if needed)..."

TMP="$EXAMPLE_DIR/example_project/input_data/genotypes/tmp"
mkdir -p "$TMP"

DATA_DIR="$EXAMPLE_DATA_DIR/genotypes"
shopt -s nullglob

for BGEN in "$DATA_DIR"/test_c*.bgen; do
  base=$(basename "$BGEN" .bgen)
  chr=$(echo "$base" | sed -E 's/.*c([0-9]+).*/\1/')
  SAMPLE="$DATA_DIR/${base}.sample"
  OUT="$TMP/c${chr}_0"

  if [[ -f "${OUT}.pgen" && ( -f "${OUT}.pvar" || -f "${OUT}.pvar.zst" ) && -f "${OUT}.psam" ]]; then
    echo "[skip] Already have $OUT"
  else
    echo "[plink2] Converting $base -> $OUT"
    plink2 --bgen "$BGEN" ref-first \
           --sample "$SAMPLE" \
           --make-pgen \
           --out "$OUT"
  fi
done

# If plink2 produced .pvar.zst, optionally decompress for readability/tools that expect plain .pvar
need_zstd=false
for z in "$TMP"/*.pvar.zst; do
  p="${z%.zst}"
  if [[ -f "$z" && ! -f "$p" ]]; then need_zstd=true; fi
done

if $need_zstd; then
  echo ""
  echo "[zstd] Decompressing .pvar.zst files..."
  if ! command -v zstd >/dev/null 2>&1 && ! command -v unzstd >/dev/null 2>&1; then
    echo "ERROR: zstd not found. Install it (brew install zstd) and rerun."
    exit 1
  fi
  for z in "$TMP"/*.pvar.zst; do
    p="${z%.zst}"
    [[ -f "$z" && ! -f "$p" ]] && zstd -d -f "$z" -o "$p"
  done
fi

echo "--- tmp contents ---"
ls -lh "$TMP" | head || true

# -------------------------------
# 3) Run pipeline (use env python)
# -------------------------------
echo ""
echo "[3/4] Running example pipeline..."

cd "$EXAMPLE_DIR"

# IMPORTANT:
# macOS often defaults to "spawn". If your code/library uses multiprocessing,
# spawn will re-import the main script and can crash unless properly guarded.
# We set start method early in a wrapper and then run the script as __main__.
conda run -n "$AGENTP_ENV" python - <<'PY'
import multiprocessing as mp
import runpy, os

try:
    mp.set_start_method("fork", force=True)
except RuntimeError:
    # start method already set; ignore
    pass

runpy.run_path("run_example.py", run_name="__main__")
PY

echo ""
echo "[4/4] Done."
echo "============================================"
echo "Example Run Complete!"
echo "============================================"
echo "Output in:"
echo "$EXAMPLE_DIR/example_project/"
