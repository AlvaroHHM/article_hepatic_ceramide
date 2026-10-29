#!/usr/bin/env bash
# Compile F1 TikZ diagram -> figures/main/F1_causal_chain.pdf
# Usage: bash analysis/figures/F1_causal_chain.sh
#
# Note: pdflatex returns non-zero on warnings even when it produces the PDF,
# so we ignore its exit code and verify the PDF exists instead.

set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
MAIN_DIR="$(cd "$SCRIPT_DIR/../../figures/main" && pwd)"
TIKZ="$MAIN_DIR/F1_causal_chain.tikz"
PDF="$MAIN_DIR/F1_causal_chain.pdf"

if [[ ! -f "$TIKZ" ]]; then
  echo "[ERROR] Source not found: $TIKZ" >&2
  exit 1
fi

cd "$MAIN_DIR"
rm -f "$PDF"
pdflatex -interaction=nonstopmode F1_causal_chain.tikz >/tmp/F1_compile.log 2>&1 || true
pdflatex -interaction=nonstopmode F1_causal_chain.tikz >>/tmp/F1_compile.log 2>&1 || true
rm -f F1_causal_chain.aux F1_causal_chain.log F1_causal_chain.out

if [[ -f "$PDF" ]] && [[ -s "$PDF" ]]; then
  SIZE=$(stat -c%s "$PDF")
  echo "[OK] F1_causal_chain.pdf (${SIZE} bytes)"
  exit 0
else
  echo "[ERROR] PDF not produced. See /tmp/F1_compile.log" >&2
  tail -15 /tmp/F1_compile.log >&2
  exit 1
fi
