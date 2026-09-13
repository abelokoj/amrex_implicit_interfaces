#!/bin/bash
# =====================================================================
# Slurm job array for the FSU RCC HPC: one array task per case.
# =====================================================================
# The sweep is a set of independent calculations, so the cases should be advanced concurrently rather than one after another. This script submits the case file as a job array, in which task N advances line N and nothing else. Where the sequential runner takes the sum of the case times, the array takes the longest single case plus whatever the tasks spend queued.
#
# Submit from amrex_implicit_interfaces/NS_fluids_lib/LSA_Rayleigh_Fortran:
#
#     scripts/submit_hpc_array.sh                  # every case in full_analysis.txt
#     scripts/submit_hpc_array.sh my_subset.txt    # a subset
#     scripts/submit_hpc_array.sh full_analysis.txt 8   # 8 ranks per case
#
# This is a wrapper rather than a submit script: it counts the cases first, because the array range has to be known at submission, and then calls sbatch with that range. Merge the per-task summaries afterwards with scripts/merge_summaries.sh.
#
# Twin pairs are safe here. Each case in the file is generated, advanced and reduced by a single task, so a perturbed calculation and its unperturbed twin are never split across tasks, machines or builds.
# ---------------------------------------------------------------------
set -euo pipefail

# Submitting from inside an existing allocation leaks that job's Slurm environment into this one, because --export=ALL carries it. A nested srun then tries to create its step in the outer job rather than in its own and is cancelled at once, which presents as tasks that end in seconds having produced nothing. Submit from a login node.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  echo "error: this appears to be running inside Slurm job ${SLURM_JOB_ID}." >&2
  echo "Submit the array from a login node instead: 'exit' this allocation, then run this script again." >&2
  echo "Set LSA_ALLOW_NESTED=1 to override, which is rarely the right choice." >&2
  [[ -n "${LSA_ALLOW_NESTED:-}" ]] || exit 1
fi


HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HERE"

CASEFILE="${1:-full_analysis.txt}"
NPROCS="${2:-8}"
ACCOUNT="${ACCOUNT:-genacc_q}"
WALLTIME="${WALLTIME:-24:00:00}"

if [[ ! -f "$CASEFILE" ]]; then
  echo "error: case file '$CASEFILE' not found in $HERE" >&2
  exit 1
fi

# Count the cases. Blank lines and comments are skipped, matching run_all.sh.
NCASE=$(grep -cve '^[[:space:]]*$' -e '^[[:space:]]*#' "$CASEFILE" || true)
if (( NCASE < 1 )); then
  echo "error: '$CASEFILE' contains no cases" >&2
  exit 1
fi

# Concurrent task limit. Submitting a hundred tasks at once is discourteous on a shared cluster and does not finish sooner, since the scheduler grants what is free either way.
THROTTLE="${THROTTLE:-8}"

# The account limits are ceilings, not defaults: a request above the limit is rejected at submission rather than trimmed. Catching that here saves a puzzling error from sbatch.
case "$ACCOUNT" in
  quicktest)          LIMIT_MIN=10 ;;
  backfill|backfill2) LIMIT_MIN=240 ;;
  *)                  LIMIT_MIN=0 ;;    # genacc_q and owner accounts, effectively unlimited here
esac
if (( LIMIT_MIN > 0 )); then
  # WALLTIME is HH:MM:SS or D-HH:MM:SS; reduce it to minutes for the comparison.
  req_min=$(awk -F'[:-]' '{ if (NF==4) print (($1*24+$2)*60)+$3; else print ($1*60)+$2 }' <<< "$WALLTIME")
  if (( req_min > LIMIT_MIN )); then
    echo "error: WALLTIME=$WALLTIME exceeds the $LIMIT_MIN minute limit of account '$ACCOUNT'." >&2
    echo "Slurm rejects such a request rather than reducing it. Lower WALLTIME, or use ACCOUNT=genacc_q, which permits up to fourteen days." >&2
    exit 1
  fi
fi

# Optional packing onto one node. By default the scheduler places the ranks wherever cores are free, which starts soonest but scatters a small job across several nodes; every rank then communicates over the network at each of the many pressure solves per step. Setting PACK=1 asks for all ranks on one node, which is faster to run and slower to start. The FSU documentation advises against constraining node counts in general, so this is off by default and worth trying only if the scattered placement proves slow.
PACK_OPT=""
if [[ -n "${PACK:-}" ]]; then
  PACK_OPT="--ntasks-per-node=$NPROCS"
  echo "  packing all $NPROCS rank(s) onto one node"
fi

mkdir -p runs/array

# Built once, here, before any task is submitted. Every task of the array shares this working directory, so a build inside a task would delete and rewrite bin/ and build/ while sibling tasks were executing those same binaries. That failure is quiet: the task exits, the sweep runner logs the case as failed, and sacct reports the task COMPLETED.
echo "building the analysis tools"
if ! make >/dev/null; then
  echo "error: build failed; nothing submitted" >&2
  exit 1
fi
for prog in make_inputs run_lsa dump_amplitude; do
  [[ -x "bin/$prog" ]] || { echo "error: bin/$prog missing after build" >&2; exit 1; }
done

echo "submitting $NCASE case(s) from $CASEFILE"
echo "  $NPROCS rank(s) per case, at most $THROTTLE at once"
echo "  account $ACCOUNT, wall-clock limit $WALLTIME per task"

sbatch --array="1-${NCASE}%${THROTTLE}" \
       --job-name="lsa_$(basename "$CASEFILE" .txt)" \
       --ntasks="$NPROCS" $PACK_OPT \
       --time="$WALLTIME" \
       --account="$ACCOUNT" \
       --output="runs/array/task-%A_%a.out" \
       --error="runs/array/task-%A_%a.err" \
       --mail-type=END,FAIL \
       --export=ALL,LSA_CASEFILE="$CASEFILE" \
       "$HERE/scripts/_array_task.sh"
