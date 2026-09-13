#!/bin/bash
# =====================================================================
# Body of one array task. Submitted by scripts/submit_hpc_array.sh; not intended to be run by hand.
# =====================================================================
# The task extracts its own line from the case file, writes it to a private one-line case file, and hands that to the ordinary sweep runner. The runner is therefore the same code path in the array as it is on a workstation, which is deliberate: a parallel submission route that reimplemented the running of a case would be a second thing to keep correct.
#
# The log and the summary are redirected to per-task files, because several tasks run at once and appending to one summary from all of them invites interleaved lines. scripts/merge_summaries.sh combines them afterwards.
# ---------------------------------------------------------------------
set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

# Any step-level variables inherited from a shell that was itself inside an allocation would misdirect the srun launches below into another job's step. Slurm sets this job's own variables, so removing the step-level ones is safe and removes an inherited value that is not.
unset SLURM_STEP_ID SLURM_STEPID SLURM_STEP_NUM_TASKS SLURM_STEP_NODELIST
unset SLURM_PMI_FD SLURM_PMI_KVS_NO_DUP_KEYS SLURM_SRUN_COMM_HOST SLURM_SRUN_COMM_PORT

CASEFILE="${LSA_CASEFILE:-full_analysis.txt}"
TASK="${SLURM_ARRAY_TASK_ID:?this script runs only as a Slurm array task}"
status=0

# Load the toolchain the solver was built with. A solver built against one MPI implementation and launched under another fails at run time rather than at link time, so these must match what 'module list' reported on the node where the solver was compiled.
# The module names are those confirmed on the login node with 'module avail'. They must match what was loaded when the solver was compiled, since a solver built against one MPI and launched under another fails at run time rather than at link time. Drop the version suffix to take the current default.
module purge
module load gnu/13
module load openmpi
module load python

mkdir -p runs/array

# The task's own line, counting only cases, so that the numbering matches what submit_hpc_array.sh counted.
LINE=$(grep -ve '^[[:space:]]*$' -e '^[[:space:]]*#' "$CASEFILE" | sed -n "${TASK}p")
if [[ -z "$LINE" ]]; then
  echo "task $TASK has no corresponding case in $CASEFILE" >&2
  exit 1
fi

MYCASE="runs/array/case_${SLURM_ARRAY_JOB_ID}_${TASK}.txt"
printf '%s\n' "$LINE" > "$MYCASE"

echo "task $TASK of array $SLURM_ARRAY_JOB_ID, $SLURM_NTASKS rank(s)"
echo "case: $LINE"
echo "started $(date)"

# Confirmed rather than assumed, since a login-node environment does not necessarily carry over to a compute node, and a failure of the Python stack costs seconds now against a whole task later. Nothing is installed here: every task shares one working directory, and concurrent installs would interfere.
if ! python3 -c 'import numpy, scipy, yt' 2>/dev/null; then
  echo "error: numpy, scipy and yt must import. Install them once before submitting; do not install inside an array task." >&2
  exit 1
fi

# DECKDIR is private to this task. Deck filenames encode the wavenumber and the amplitude but not the resolution or the domain size, so the resolution study and the boundary case would otherwise all write inputs.k0.7.eps0.02 into one directory and overwrite one another, and a task could advance a deck another task had just replaced.
LAUNCHER=srun \
DECKDIR="runs/array/decks_${SLURM_ARRAY_JOB_ID}_${TASK}" \
LOG="runs/array/run_${SLURM_ARRAY_JOB_ID}_${TASK}.log" \
SUMMARY="runs/array/summary_${SLURM_ARRAY_JOB_ID}_${TASK}.tsv" \
  ./scripts/run_all.sh "$MYCASE" "$SLURM_NTASKS" || status=$?

echo "finished $(date), status $status"
# Propagated so that sacct reports a failed task as failed. Without this the task exits 0 and a sweep in which every case died is recorded as COMPLETED.
exit $status
