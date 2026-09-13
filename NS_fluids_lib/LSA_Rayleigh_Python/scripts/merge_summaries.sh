#!/bin/bash
# =====================================================================
# Merge the per-task summaries written by a Slurm array run.
# =====================================================================
# Each array task writes its own summary so that concurrent tasks do not interleave lines in a shared file. This combines them into runs/summary.tsv, in the same format the sequential runner produces, so that the plotting and reporting code need not know which route produced the results.
#
# Run from amrex_implicit_interfaces/NS_fluids_lib/LSA_Rayleigh_Python:
#
#     ./scripts/merge_summaries.sh
# ---------------------------------------------------------------------
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HERE"

SRC="runs/array"
OUT="runs/summary.tsv"

shopt -s nullglob
parts=("$SRC"/summary_*.tsv)
if (( ${#parts[@]} == 0 )); then
  echo "no per-task summaries under $SRC; nothing to merge" >&2
  exit 1
fi

# Header from the first part, then every data line from all of them, sorted by label so that the merged file is reproducible regardless of the order the tasks happened to finish in.
head -n 1 "${parts[0]}" > "$OUT"
tail -q -n +2 "${parts[@]}" | sort -u >> "$OUT"

n=$(( $(wc -l < "$OUT") - 1 ))
echo "merged ${#parts[@]} task summaries into $OUT ($n case(s))"

if command -v column >/dev/null 2>&1; then
  column -t -s $'\t' "$OUT"
fi
