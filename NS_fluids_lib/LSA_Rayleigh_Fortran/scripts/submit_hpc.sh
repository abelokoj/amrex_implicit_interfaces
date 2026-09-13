#!/bin/bash
# =====================================================================
# Slurm submit script for the FSU Research Computing Center HPC.
# =====================================================================
# Submits the whole sweep as one batch job, so that every case in full_analysis.txt is advanced without further intervention.
#
# Submit from amrex_implicit_interfaces/NS_fluids_lib/LSA_Rayleigh_Fortran:
#
#     sbatch scripts/submit_hpc.sh
#     sbatch scripts/submit_hpc.sh my_subset.txt
#
# Then: squeue --me to watch it, scancel JOBID to stop it, and tail runs/run_all.log to follow progress.
#
# An interactive job through srun --pty is the wrong instrument for this. It ties the calculation to a live terminal, so a dropped connection ends it, and the interactive queues carry short time limits. Interactive work is for the smoke test of Step 4 and for inspecting output, not for a sweep that runs for hours or days.
#
# Two choices below deserve comment. The launcher is srun rather than mpirun, which is what the RCC documentation recommends, because srun inherits the allocation from the scheduler; this also avoids the singleton failure mode in which a launcher that cannot reach a process manager starts independent single-rank copies instead of one job. And the number of ranks is taken from SLURM_NTASKS rather than written twice, so that the value requested from the scheduler and the value passed to the solver cannot drift apart.
# ---------------------------------------------------------------------
#SBATCH --job-name=lsa_rayleigh
#SBATCH --ntasks=16
#SBATCH --time=48:00:00
#SBATCH --output=runs/slurm-%j.out
#SBATCH --error=runs/slurm-%j.err
#SBATCH --mail-type=END,FAIL

# Slurm account. genacc_q is the general access account and allows up to fourteen days; backfill and backfill2 are limited to four hours and backfill2 may be pre-empted, so neither suits a long sweep. Replace this with an owner-based account if the group holds one.
#SBATCH --account=genacc_q

set -euo pipefail

# Submitting from inside an existing allocation leaks that job's Slurm environment into this one, because --export=ALL carries it. A nested srun then tries to create its step in the outer job rather than in its own and is cancelled at once, which presents as tasks that end in seconds having produced nothing. Submit from a login node.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  echo "error: this appears to be running inside Slurm job ${SLURM_JOB_ID}." >&2
  echo "Submit the array from a login node instead: 'exit' this allocation, then run this script again." >&2
  echo "Set LSA_ALLOW_NESTED=1 to override, which is rarely the right choice." >&2
  [[ -n "${LSA_ALLOW_NESTED:-}" ]] || exit 1
fi


CASEFILE="${1:-full_analysis.txt}"

# The job begins in the directory from which it was submitted, which the instructions above fix as LSA_Rayleigh_Fortran. Every path below is relative to it, and every file the job writes stays inside it.
mkdir -p runs

# Load the same toolchain the solver was built with. A solver built against one MPI implementation and launched under another is a common cause of failures that appear at run time rather than at link time, so these must match what was used for 'make -j16' in NS_fluids_lib. Adjust the module names to whatever 'module list' reports on the build node.
# The module names are those confirmed on the login node with 'module avail'. They must match what was loaded when the solver was compiled, since a solver built against one MPI and launched under another fails at run time rather than at link time. Drop the version suffix to take the current default.
module purge
module load gnu/13
module load openmpi

echo "job $SLURM_JOB_ID on $SLURM_JOB_NUM_NODES node(s), $SLURM_NTASKS task(s)"
echo "started $(date)"
echo "case file: $CASEFILE"

# Build the analysis tools on the compute node rather than assuming a build from the login node is present and current.
make

# LAUNCHER is read by run_all.sh in place of mpirun. SLURM_NTASKS carries the rank count granted by the scheduler.
LAUNCHER=srun ./scripts/run_all.sh "$CASEFILE" "$SLURM_NTASKS"

echo "finished $(date)"
echo "summary written to runs/summary.tsv"
