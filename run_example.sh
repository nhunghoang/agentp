#!/bin/bash
# run_example.sh - Auto-detecting example runner for dev_agentp
# Works on macOS ARM, macOS Intel, Linux
# Usage: bash run_example.sh

set -euo pipefail

echo "============================================"
echo "Running dev_agentp Example"
echo "============================================"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$SCRIPT_DIR"
EX_TAR="$SCRIPT_DIR/example.tar.gz"
BOX_URL="https://vanderbilt.box.com/shared/static/o5s5bwwey0clt8ktabqg1abmubg2t7in.gz"

# ------------------------------------------------------------
# 1. AUTO-DETECT BGENIX
# ------------------------------------------------------------
echo ""
echo "[Detect] Searching for bgenix..."

BGENIX_CANDIDATES=(
  "$HOME/miniconda3/envs/bgen/bin/bgenix"
  "$HOME/miniforge3/envs/bgen/bin/bgenix"
  "$HOME/miniconda3/envs/bgenix/bin/bgenix"
  "$HOME/miniforge3/envs/bgenix/bin/bgenix"
  "$HOME/micromamba/envs/bgen/bin/bgenix"
  "$HOME/.local/bin/bgenix"
)

BGENIX_BIN=""

for p in "${BGENIX_CANDIDATES[@]}"; do
  if [[ -x "$p" ]]; then
    BGENIX_BIN="$p"
    break
  fi
done

# If missing, try PATH
if [[ -z "$BGENIX_BIN" ]] && command -v bgenix &>/dev/null; then
  BGENIX_BIN="$(command -v bgenix)"
fi

if [[ -z "$BGENIX_BIN" ]]; then
  echo "ERROR: bgenix not found. Please run install.sh."
  exit 1
fi

export PATH="$(dirname "$BGENIX_BIN"):$PATH"
echo "[OK] Using bgenix: $BGENIX_BIN"

# ------------------------------------------------------------
# 2. AUTO-DETECT PLINK2
# ------------------------------------------------------------
echo ""
echo "[Detect] Searching for plink2..."

if command -v plink2 &>/dev/null; then
  PLINK2_BIN="$(command -v plink2)"
elif [[ -x "$HOME/.local/bin/plink2" ]]; then
  PLINK2_BIN="$HOME/.local/bin/plink2"
  export PATH="$HOME/.local/bin:$PATH"
else
  echo "ERROR: plink2 not found. Please run install.sh."
  exit 1
fi

echo "[OK] Using plink2: $PLINK2_BIN"

# ------------------------------------------------------------
# 3. DOWNLOAD & EXTRACT EXAMPLE DATA
# ------------------------------------------------------------
echo ""
echo "[1/4] Downloading example data..."

if command -v curl &>/dev/null; then
  curl -fL "$BOX_URL" -o "$EX_TAR"
elif command -v wget &>/dev/null; then
  wget -O "$EX_TAR" "$BOX_URL"
else
  echo "ERROR: curl or wget required."
  exit 1
fi

echo "Verifying archive..."
tar -tzf "$EX_TAR" >/dev/null

echo "Extracting..."
tar -xzf "$EX_TAR" -C "$SCRIPT_DIR"
rm -f "$EX_TAR"

# ------------------------------------------------------------
# 4. LINK example_data
# ------------------------------------------------------------
echo ""
echo "[2/4] Setting up example data..."

FOUND="$(find "$SCRIPT_DIR" -maxdepth 3 -type d -name 'example_data' -not -path "$REPO/example/example_data" | head -n1 || true)"
[[ -z "$FOUND" ]] && { echo "ERROR: example_data not found."; exit 1; }

rm -rf "$REPO/example/example_data"
ln -s "$FOUND" "$REPO/example/example_data"

echo "example_data -> $(readlink "$REPO/example/example_data")"

# ------------------------------------------------------------
# 5. CONVERT BGEN → PGEN
# ------------------------------------------------------------
echo ""
echo "[3/4] Converting genotype files..."

TMP="$REPO/example/example_project/input_data/genotypes/tmp"
mkdir -p "$TMP"

DATA_DIR="$REPO/example/example_data/genotypes"
shopt -s nullglob

for BGEN in "$DATA_DIR"/test_c*.bgen; do
  base=$(basename "$BGEN" .bgen)
  chr=$(echo "$base" | sed -E 's/.*c([0-9]+).*/\1/')
  SAMPLE="$DATA_DIR/${base}.sample"
  OUT="$TMP/c${chr}_0"

  if [[ -f "${OUT}.pgen" && ( -f "${OUT}.pvar" || -f "${OUT}.pvar.zst" ) && -f "${OUT}.psam" ]]; then
    echo "[skip] Already have $OUT"
  else
    echo "[plink2] Converting $base → $OUT"
    plink2 --bgen "$BGEN" ref-first \
           --sample "$SAMPLE" \
           --make-pgen \
           --out "$OUT"
  fi
done

# ------------------------------------------------------------
# 6. HANDLE zstd DECOMPRESSION (if needed)
# ------------------------------------------------------------
need_zstd=false
for z in "$TMP"/*.pvar.zst; do
  p="${z%.zst}"
  if [[ -f "$z" && ! -f "$p" ]]; then need_zstd=true; fi
done

if $need_zstd; then
  echo ""
  echo "[zstd] Decompressing .pvar.zst files..."

  if ! command -v zstd &>/dev/null && ! command -v unzstd &>/dev/null; then
    if [[ "$OSTYPE" == "darwin"* ]]; then
      brew install zstd
    else
      sudo apt-get update -qq
      sudo apt-get install -y zstd
    fi
  fi

  for z in "$TMP"/*.pvar.zst; do
    p="${z%.zst}"
    [[ -f "$z" && ! -f "$p" ]] && zstd -d -f "$z" -o "$p"
  done
fi

echo "--- tmp contents ---"
ls -lh "$TMP" | head || true

# ------------------------------------------------------------
# 7. AUTO-DETECT PYTHON WITH agentp INSTALLED
# ------------------------------------------------------------
echo ""
echo "[Detect] Searching for Python with agentp installed..."

PY_CANDIDATES=(
  "$HOME/miniconda3/bin/python"
  "$HOME/miniforge3/bin/python"
  "$HOME/micromamba/bin/python"
  "$(command -v python3 || true)"
  "$(command -v python || true)"
)

AGENTP_PY=""

for py in "${PY_CANDIDATES[@]}"; do
  if [[ -x "$py" ]]; then
    if "$py" - << 'PY' 2>/dev/null
import agentp
print("FOUND")
PY
    then
      AGENTP_PY="$py"
      break
    fi
  fi
done

if [[ -z "$AGENTP_PY" ]]; then
  echo "ERROR: Could not find Python installation containing agentp."
  echo "Please run install.sh."
  exit 1
fi

echo "[OK] Using Python: $AGENTP_PY"

# ------------------------------------------------------------
# 8. RUN THE PYTHON PIPELINE
# ------------------------------------------------------------
echo ""
echo "[4/4] Running example pipeline..."
cd "$REPO/example"

"$AGENTP_PY" run_example.py || {
  echo "WARNING: Pipeline returned an error. Check logs above."
}

echo ""
echo "============================================"
echo "Example Run Complete!"
echo "============================================"
echo "Output in:"
echo "$REPO/example/example_project/"
