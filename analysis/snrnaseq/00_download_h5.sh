#!/bin/bash
# ============================================================
# 00_download_h5.sh
# Download the 28 H5 files from GSE250589 (Zheng et al. 2025).
# Required by snrnaseq scripts 05-10 (MOL3_EAE analyses).
#
# Size: ~2.3 GB. Time: 10-30 min depending on connection.
# Usage: bash analysis/snrnaseq/00_download_h5.sh
# ============================================================

set -eu

# Move to repo root
cd "$(dirname "$0")/../.."

H5_DIR="data/GSE250589/h5"
URL_RAW="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250589/suppl/GSE250589_RAW.tar"
RAW_TAR="/tmp/GSE250589_RAW.tar"

mkdir -p "$H5_DIR"

# 1. Check if already downloaded
N_H5=$(find "$H5_DIR" -maxdepth 1 -name "*.h5" 2>/dev/null | wc -l | tr -d ' ')
if [ "$N_H5" -ge 28 ]; then
    echo "[SKIP] $N_H5 H5 files already in $H5_DIR"
    exit 0
fi
echo "[INFO] Found $N_H5/28 H5. Downloading the rest (~2.3 GB)."

# 2. Download the raw tar
if [ ! -f "$RAW_TAR" ]; then
    echo "[INFO] Downloading GSE250589_RAW.tar..."
    if command -v wget >/dev/null 2>&1; then
        wget -O "$RAW_TAR" "$URL_RAW"
    elif command -v curl >/dev/null 2>&1; then
        curl -L -o "$RAW_TAR" "$URL_RAW"
    else
        echo "[ERROR] Need wget or curl" >&2
        exit 1
    fi
fi

# 3. Extract H5 files
TMP_EXTRACT=$(mktemp -d)
trap "rm -rf $TMP_EXTRACT" EXIT
echo "[INFO] Extracting H5 files to $H5_DIR ..."
tar -xf "$RAW_TAR" -C "$TMP_EXTRACT"
find "$TMP_EXTRACT" -maxdepth 2 -name "*.h5" -exec mv {} "$H5_DIR/" \;

# 4. Verify
N_H5=$(find "$H5_DIR" -maxdepth 1 -name "*.h5" 2>/dev/null | wc -l | tr -d ' ')
if [ "$N_H5" -lt 28 ]; then
    echo "[WARN] Only $N_H5 H5 extracted (expected 28)" >&2
    exit 1
fi
echo "[OK] $N_H5 H5 files ready in $H5_DIR"
