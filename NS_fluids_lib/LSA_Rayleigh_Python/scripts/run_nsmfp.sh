#!/usr/bin/env bash
# =====================================================================
# Rayleigh-Plateau NS-MFP linear stability analysis -- driver script
#                        (Python implementation)
# =====================================================================
# Advances the unperturbed twin and the perturbed case for one wavenumber,
# then post-processes the pair into a growth rate.
#
# Run from amrex_implicit_interfaces/NS_fluids_lib/LSA_Rayleigh_Python:
#
#     ./scripts/run_nsmfp.sh [kr0] [eps] [nprocs] [cells_per_r0] [solver]
#
# All paths are resolved relative to the solver tree, so nothing outside it
# need be supplied. The solver defaults to
#
#     amrex_implicit_interfaces/NS_fluids_lib/amr2d.gnu.MPI.ex
#
# which is where 'make -j16' in NS_fluids_lib places it. Every file this script writes stays inside this directory, under ./runs, so that the working tree of the solver is left untouched.
#
# Example:
#     ./scripts/run_nsmfp.sh 0.7 2e-2 8            # the reference case, 8 cells per radius
#     ./scripts/run_nsmfp.sh 0.7 2e-2 8 16         # the same case at 16 cells per radius
#
# The fourth argument is the radial resolution, and it dominates the cost. Work scales roughly as the cube of it, because the capillary stability limit on the time step varies as dx^(3/2) while the cell count varies as dx^(-2). The default of 8 is the reference resolution of the study and completes in minutes; 16 costs approximately eight times as much and 32 approximately sixty times, so neither should be started before the reference case has been seen to work.
#
# CRITICAL: the twin and the perturbed case are launched with the SAME
# number of MPI ranks. The twin subtraction stands in for the NS-MFP body
# force B_f and presumes an identical discretisation; altering the domain
# decomposition between them perturbs round-off-level behaviour, which is
# the scale at which the measurement is made.
# ---------------------------------------------------------------------
set -euo pipefail

KR0="${1:-0.7}"
EPS="${2:-2e-2}"
NPROCS="${3:-4}"
CELLS="${4:-8}"

# Launcher. Set LAUNCHER=srun under Slurm, where srun inherits the allocation and integrates with the scheduler.
LAUNCHER="${LAUNCHER:-mpirun}"
if [[ "$LAUNCHER" == *srun* ]]; then
  # See scripts/run_all.sh for why srun needs its binding disabled here.
  LAUNCHER_OPTS="${LAUNCHER_OPTS:---cpu-bind=none}"
else
  LAUNCHER_OPTS="${LAUNCHER_OPTS:-}"
fi

# ROOT is LSA_Rayleigh_Python; NSDIR is NS_fluids_lib, which holds both the
# solver executable and the run directories.
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NSDIR="$(cd "$ROOT/.." && pwd)"
PY="$ROOT/python"
RUNS="$ROOT/runs"
DECKS="$RUNS/decks"

EXE="${5:-$NSDIR/amr2d.gnu.MPI.ex}"

if [[ ! -x "$EXE" ]]; then
  echo "error: solver '$EXE' is not executable." >&2
  echo "Build it with 'cd $NSDIR && make -j16'; see README.md." >&2
  exit 1
fi
EXE="$(cd "$(dirname "$EXE")" && pwd)/$(basename "$EXE")"

if ! python3 -c 'import numpy, scipy' 2>/dev/null; then
  echo "error: numpy and scipy are required; see requirements.txt." >&2
  exit 1
fi


# A launcher that cannot reach its process manager silently starts N independent
# single-rank copies instead of one N-rank job. Each copy would then advance the
# whole problem in the same directory and overwrite the others' output, which is
# far worse than a crash because it looks like success. Detect that here.
if (( NPROCS > 1 )); then
  probe=0
  if [[ "$LAUNCHER" == mpirun* ]]; then
    probe=$(mpirun -n 2 hostname 2>&1 | grep -c 'singletons will be started' || true)
  fi
  # A launcher from a different MPI than the solver was linked against fails outright rather than degrading, and the message names a PMI version rather than anything obviously to do with the launcher.
  if ! $LAUNCHER -n 1 "$EXE" --version >/dev/null 2>&1; then
    if $LAUNCHER -n 1 "$EXE" 2>&1 | grep -q 'unsupported PMI version'; then
      echo "error: '$LAUNCHER' comes from a different MPI installation than the solver was linked against." >&2
      echo "Compare 'mpirun --version' with 'ldd $EXE | grep -i mpi' and use the launcher belonging to that library; see RUNNING guide, Section 15." >&2
      exit 1
    fi
  fi
  if (( probe > 0 )); then
    echo "error: mpirun cannot reach a process manager on this machine and would start $NPROCS independent single-rank copies of the solver, all writing into the same directory." >&2
    echo "Re-run with 1 rank, or repair the launcher; see Running_the_Analysis_Python.md, Section 11." >&2
    exit 1
  fi
fi

echo "=== 1/4  generating matched decks (k r0=$KR0, eps=$EPS, $CELLS cells per radius) ==="
# Run from the project root so that the default template path resolves.
( cd "$ROOT" && python3 "$PY/make_inputs.py" single --kr0 "$KR0" \
    --eps "$EPS" --cells-per-r0 "$CELLS" --rmax 4.0 --n-periods 5.0 \
    --snapshots 120 --outdir "$DECKS" )

# make_inputs.py names decks with the '%g' form of each value, e.g. eps 2e-2
# becomes "0.02"; the tags are computed the same way here so that they agree.
EPS_TAG="$(python3 -c 'import sys; print(f"{float(sys.argv[1]):g}")' "$EPS")"
KR0_TAG="$(python3 -c 'import sys; print(f"{float(sys.argv[1]):g}")' "$KR0")"

BASE_DECK="$DECKS/inputs.k${KR0_TAG}.base"
PERT_DECK="$DECKS/inputs.k${KR0_TAG}.eps${EPS_TAG}"
BASE_DIR="$RUNS/k${KR0_TAG}_base"
PERT_DIR="$RUNS/k${KR0_TAG}_eps${EPS_TAG}"

for d in "$BASE_DECK" "$PERT_DECK"; do
  [[ -f "$d" ]] || { echo "error: expected deck '$d' was not generated" >&2; exit 1; }
done

run_case () {
  local deck="$1" dir="$2" label="$3"
  echo
  echo "--- running $label -> $dir"
  # Each case runs in its own directory because the solver writes its
  # plotfiles and checkpoints into the current working directory.
  rm -rf "$dir"; mkdir -p "$dir"
  ( cd "$dir" && cp "$deck" ./inputs \
      && $LAUNCHER $LAUNCHER_OPTS -n "$NPROCS" "$EXE" inputs \
           > run.out 2> run.err )
  local n
  # The solver does not honour amr.plot_file: with ns.visual_nddata_format=1
  # the plotfile directories are named nddataPLT<step>, not plt<step>.
  n=$(find "$dir" -maxdepth 1 -type d -name 'nddataPLT*' | wc -l)
  echo "    wrote $n plotfiles"
  if (( n < 2 )); then
    echo "error: $label produced fewer than 2 plotfiles." >&2
    echo "Inspect $dir/run.err and confirm ns.visual_nddata_format=1." >&2
    exit 1
  fi
}

echo
echo "=== 2/4  running the twin pair on $NPROCS rank(s) ==="
run_case "$BASE_DECK" "$BASE_DIR" "unperturbed twin (base state)"
run_case "$PERT_DECK" "$PERT_DIR" "perturbed (eps=$EPS)"

echo
echo "=== 3/4  extracting the growth rate ==="
python3 "$PY/run_lsa.py" analyse \
  --perturbed "$PERT_DIR" \
  --base "$BASE_DIR" \
  --kr0 "$KR0" \
  --skip 60 \
  --save "$RUNS/result_k${KR0_TAG}_eps${EPS_TAG}.json"

echo
echo "=== 4/4  exporting the figure data ==="
# Exported now, while the plotfiles are present, so that all four figures of this case can be redrawn later from a few hundred kilobytes of plain text without yt and without the plotfiles. This is what makes the figures reproducible on a workstation when the calculation ran on a cluster.
BUNDLE="$RUNS/bundles/k${KR0_TAG}_eps${EPS_TAG}"
mkdir -p "$BUNDLE"
LAMBDA=$(awk '/^ns.yblob/{print $3; exit}' "$PERT_DIR/inputs")
python3 "$PY/dump_amplitude.py" "$PERT_DIR" "$LAMBDA" "$BUNDLE/amplitude.dat"
python3 "$PY/dump_profiles.py" "$PERT_DIR" "$LAMBDA" "$BUNDLE" --field y_velocity

echo
echo "Draw the figures with:"
echo "  python3 python/plot_results.py --bundle $BUNDLE --kr0 $KR0"
echo
echo "Done. Before this number is quoted, run the proportionality test:"
echo "  see Running_the_Analysis_Python.md, Step 7."
