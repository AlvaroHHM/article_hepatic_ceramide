#!/bin/bash
# ============================================================
# 00_download_h5.sh
# Download the 28 H5 files from GSE250589 (Zheng et al. 2025).
# Required by snrnaseq scripts 05-10 (MOL3_EAE analyses).
#
# Uses STREAMING extraction: downloads the GEO RAW.tar via a
# pipe and extracts only *.h5 files, so the full ~40 GB tar
# never touches the disk. Only ~2.3 GB of H5 files remain.
#
# Size: ~2.3 GB extracted (no 40 GB temp).
# Time: 10-30 min depending on connection.
# Usage: bash analysis/snrnaseq/00_download_h5.sh
# ============================================================

set -eu
cd "$(dirname "$0")/../.."

H5_DIR="data/GSE250589/h5"
URL_RAW="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250589/suppl/GSE250589_RAW.tar"

mkdir -p "$H5_DIR"

N_H5=$(find "$H5_DIR" -maxdepth 1 -name "*.h5" 2>/dev/null | wc -l | tr -d ' ')
if [ "$N_H5" -ge 28 ]; then
    echo "[SKIP] $N_H5 H5 files already in $H5_DIR"
    exit 0
fi

echo "[INFO] Found $N_H5/28 H5. Streaming extraction from GEO FTP."
echo "[INFO] This extracts only *.h5 (~2.3 GB), not the full tar (~40 GB)."

# Stream the tar through a pipe and extract only *.h5
if command -v wget >/dev/null 2>&1; then
    wget -O - "$URL_RAW" 2>/dev/null | tar -xf - -C "$H5_DIR" --wildcards '*.h5'
elif command -v curl >/dev/null 2>&1; then
    curl -sL "$URL_RAW" | tar -xf - -C "$H5_DIR" --wildcards '*.h5'
else
    echo "[ERROR] Need wget or curl" >&2
    exit 1
fi

N_H5=$(find "$H5_DIR" -maxdepth 1 -name "*.h5" 2>/dev/null | wc -l | tr -d ' ')
if [ "$N_H5" -lt 28 ]; then
    echo "[WARN] Only $N_H5 H5 extracted (expected 28)" >&2
    echo "[WARN] Check $H5_DIR manually" >&2
    exit 1
fi
echo "[OK] $N_H5 H5 files ready in $H5_DIR"
