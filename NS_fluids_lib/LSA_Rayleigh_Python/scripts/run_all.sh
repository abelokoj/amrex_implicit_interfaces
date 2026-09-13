#!/usr/bin/env bash
# =====================================================================
# run_all.sh -- advance every calculation listed in a sweep file, then extract
# the interface mode amplitude from each.
#
# Run from amrex_implicit_interfaces/NS_fluids_lib/LSA_Rayleigh_Python:
#
#     ./scripts/run_all.sh                          # full_analysis.txt, 1 rank
#     ./scripts/run_all.sh full_analysis.txt 16   # 16 ranks
#     ./scripts/run_all.sh my_subset.txt 8 ../amr2d.gnu.MPI.ex
#

# The launcher is selected by the LAUNCHER environment variable, which defaults to mpirun. Under Slurm it should be set to srun, which the FSU Research Computing Center documentation recommends in preference to mpirun because srun integrates with the scheduler and inherits the allocation, so the rank count and the placement come from the job rather than from a command-line argument that may disagree with it. Setting it also avoids the singleton failure mode, in which a launcher that cannot reach a process manager starts independent single-rank copies instead of one job:
#
#     LAUNCHER=srun ./scripts/run_all.sh full_analysis.txt 16
#
# The script is designed to be started once and left alone. It is also safe to
# interrupt and restart: a case whose amplitude record already exists is
# skipped, and a case that is partly advanced resumes from its most recent
# checkpoint. Progress is appended to runs/run_all.log, and a machine
# readable summary is written to runs/summary.tsv as each case finishes.
#
# Cases are attempted in file order and a failure does not stop the sweep, so
# that one unstable configuration cannot cost the whole run. Failures are
# recorded in the log and marked in the summary.
# =====================================================================
set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HERE"

CASEFILE="${1:-full_analysis.txt}"
NPROCS="${2:-1}"
EXE="${3:-../amr2d.gnu.MPI.ex}"

# Launcher. Set LAUNCHER=srun under Slurm; see the header of this script.
LAUNCHER="${LAUNCHER:-mpirun}"
if [[ "$LAUNCHER" == mpirun* ]]; then
  LAUNCHER_OPTS="${LAUNCHER_OPTS:---allow-run-as-root --oversubscribe}"
else
  # --cpu-bind=none is required, not cosmetic. Slurm satisfies a request for N tasks from whatever cores are free, which on a busy cluster means a few cores on each of several nodes, in fragmented masks such as 0x140000 or 0x000001. srun's default binding cannot be mapped onto those and every launch is refused with "Unable to satisfy cpu bind request", so the solver never starts. Leaving placement to the operating system costs a little locality and avoids the failure entirely.
  LAUNCHER_OPTS="${LAUNCHER_OPTS:---cpu-bind=none}"
fi

RUNROOT="runs"
# The deck directory may be redirected. Generated deck names encode the wavenumber and amplitude but not the resolution or the domain size, so two cases that differ only in those would be written to the same filename; concurrent jobs must therefore each be given their own directory. See scripts/_array_task.sh.
DECKDIR="${DECKDIR:-$RUNROOT/decks}"
RECDIR="$RUNROOT/records"
# Figure data. Each case gets its own bundle of interface profiles and field slices, written straight after its amplitude record, so that every figure of the study can be redrawn later from plain text. The plotfiles usually stay on the machine that ran the sweep and amount to tens of gigabytes; a bundle is a few hundred kilobytes and travels in one copy. Set FIELDS empty to skip the field slices, which are the larger part.
BUNDLEDIR="${BUNDLEDIR:-$RUNROOT/bundles}"
FIELDS="${FIELDS:-y_velocity}"
# The log and the summary may be redirected so that concurrent jobs, such as the tasks of a Slurm array, each write their own and are merged afterwards rather than appending to one file at once.
LOG="${LOG:-$RUNROOT/run_all.log}"
SUMMARY="${SUMMARY:-$RUNROOT/summary.tsv}"

# Solver tolerances. These must sit well below the smallest perturbation
# amplitude in the sweep, or the smallest case measures solver residual rather
# than physics.
# Overridable, in the manner of CHECK_INT, so that a tolerance or a stop time may be varied for one sweep without editing this file.
SOLVER_OPTS="${SOLVER_OPTS:-mac.mac_abs_tol=1.0e-10 mac.visc_abs_tol=1.0e-10}"

# Checkpoint interval. This must be smaller than the number of steps a single
# scheduler allocation can complete, or no checkpoint is ever written and each
# restart silently repeats from the beginning.
CHECK_INT="${CHECK_INT:-50}"

mkdir -p "$RUNROOT" "$DECKDIR" "$RECDIR" "$BUNDLEDIR"

log () { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

if [[ ! -f "$CASEFILE" ]]; then
  log "ERROR: case file '$CASEFILE' not found"
  exit 1
fi
if [[ ! -x "$EXE" ]]; then
  log "ERROR: solver '$EXE' is not executable; build it with"
  log "       cd $(cd .. && pwd) && make -j16"
  exit 1
fi
# Resolved to an absolute path because each case is advanced from its own run
# directory, where a path relative to this one would no longer refer to it.
EXE="$(cd "$(dirname "$EXE")" && pwd)/$(basename "$EXE")"

if ! python3 -c 'import numpy, scipy, yt' 2>/dev/null; then
  log "ERROR: numpy, scipy and yt are required; see requirements.txt"
  exit 1
fi

[[ -f "$SUMMARY" ]] || printf 'label\tkr0\teps\tcells\trmax\tstatus\tsteps\tt_final\trecord\n' > "$SUMMARY"

if (( NPROCS > 1 )); then
  if [[ "$LAUNCHER" == mpirun* ]] && \
     mpirun -n 2 hostname 2>&1 | grep -q 'singletons will be started'; then
    log "ERROR: mpirun cannot reach a process manager and would start $NPROCS independent single-rank copies of the solver, all writing into the same directory. Re-run with 1 rank or repair the launcher; see Running_the_Analysis_Python.md, Section 11."
    exit 1
  fi
fi

log "sweep '$CASEFILE' on $NPROCS rank(s), solver '$EXE'"

# ---------------------------------------------------------------------
# Advance one case to completion, restarting as often as necessary.
#
# The solver is invoked repeatedly rather than once because a single
# invocation may be cut short by a scheduler limit; each invocation resumes
# from the latest checkpoint. The loop terminates when the step count stops
# advancing, which covers both normal completion and a stalled calculation.
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
# Domain decomposition. AMReX aborts with "cannot have a dormant/idle processor on any level" when a rank is given no box, so the number of boxes must be at least the number of ranks. The box count follows from amr.n_cell and amr.max_grid_size, both of which are known once the deck exists, so the largest max_grid_size that supplies enough boxes is chosen here and passed on the command line. Where even the smallest sensible value cannot supply enough, the rank count is reduced instead, since running a case on fewer ranks is always preferable to not running it.
#
# Larger boxes are preferred because each carries a ghost region whose cost grows as the surface-to-volume ratio worsens, so the search runs from 32 downwards rather than simply taking the smallest.
# ---------------------------------------------------------------------
decompose () {
  local deck="$1" want="$2" nx nz mgs boxes cgmin smallest smallest_boxes
  nx=$(awk '/^amr.n_cell/{print $3; exit}' "$deck")
  nz=$(awk '/^amr.n_cell/{print $4; exit}' "$deck")
  if [[ -z "$nx" || -z "$nz" ]]; then echo "32 $want"; return; fi

  # The solver requires amr.max_grid_size >= cg.min_max_grid_size and aborts otherwise, so the floor is read from the deck rather than assumed. Lowering cg.min_max_grid_size to admit smaller boxes was rejected: it alters the bottom of the multigrid solve, and every case of the sweep should be advanced under the same solver settings so that the measured rates remain comparable. Reducing the rank count instead changes only the decomposition.
  cgmin=$(awk '/^cg.min_max_grid_size/{print $3; exit}' "$deck")
  [[ -z "$cgmin" ]] && cgmin=16

  smallest=""; smallest_boxes=0
  for mgs in 32 16 8; do
    (( mgs < cgmin )) && continue
    boxes=$(( ( (nx + mgs - 1) / mgs ) * ( (nz + mgs - 1) / mgs ) ))
    smallest=$mgs; smallest_boxes=$boxes
    if (( boxes >= want )); then echo "$mgs $want"; return; fi
  done

  # Not enough boxes at any permitted size: run on as many ranks as there are boxes.
  [[ -z "$smallest" ]] && { echo "$cgmin 1"; return; }
  (( smallest_boxes < 1 )) && smallest_boxes=1
  echo "$smallest $smallest_boxes"
}

advance () {
  local dir="$1" deck="$2" label="$3"
  mkdir -p "$dir"
  [[ -f "$dir/inputs" ]] || cp "$deck" "$dir/inputs"

  local mgs ranks
  read -r mgs ranks <<< "$(decompose "$dir/inputs" "$NPROCS")"
  if (( ranks < NPROCS )); then
    log "  $label: only $ranks box(es) available, running on $ranks rank(s) rather than $NPROCS" >&2
  fi

  local prev=-1 now=0 attempt=0
  while (( attempt < 200 )); do
    attempt=$(( attempt + 1 ))
    local last restart=""
    last=$(cd "$dir" && ls -d chk* 2>/dev/null | sort | tail -1)
    [[ -n "$last" ]] && restart="amr.restart=$last"

    # stdin is redirected because mpirun would otherwise consume the case
    # file that the main loop is reading, so that only the first case runs.
    ( cd "$dir" && $LAUNCHER $LAUNCHER_OPTS -n "$ranks" \
        "$EXE" inputs $SOLVER_OPTS amr.max_grid_size=$mgs \
        amr.check_int=$CHECK_INT $restart \
        < /dev/null >> solver.out 2>> solver.err )

    now=$(grep -h '^STEP' "$dir/solver.out" 2>/dev/null | tail -1 \
          | awk '{print $3}')
    now=${now:-0}
    if [[ "$now" == "$prev" ]]; then
      # Zero steps on the first attempt means the solver never started rather than that it finished. Reporting the reason here saves reading two log files to discover that the launcher, and not the analysis, is at fault.
      if (( attempt <= 2 )) && [[ "$now" == "0" ]]; then
        log "  $label: the solver produced no time steps; see $dir/solver.err" >&2
        if [[ -s "$dir/solver.err" ]]; then
          # Skipped: the X11 and launcher chatter that heads the file and says nothing about the failure.
          log "  $label: $(grep -vE 'No protocol specified|^ *$' "$dir/solver.err" \
                 | grep -iE 'abort|error|fail|expecting|not found' \
                 | head -3 | tr '\n' ' ')" >&2
        fi
      fi
      break                     # no progress: finished, or stalled
    fi
    prev="$now"
    # To stderr: this function's stdout is captured by the caller.
    log "  $label: step $now" >&2
  done
  echo "$now"
}

# ---------------------------------------------------------------------
# Main loop over the case file.
# ---------------------------------------------------------------------
while read -r label kr0 eps cells rmax nper snaps rest; do
  # Skip blank lines and comments. The test is on the first field, so an
  # indented comment is handled as well as one in the first column.
  [[ -z "${label:-}" ]] && continue
  [[ "${label:0:1}" == "#" ]] && continue

  record="$RECDIR/${label}.dat"
  if [[ -s "$record" ]]; then
    log "$label: already complete, skipping"
    continue
  fi

  log "$label: kr0=$kr0 eps=$eps cells=$cells rmax=$rmax periods=$nper"

  # A label ending in .base denotes the unperturbed twin, for which the
  # amplitude field is meaningless and is forced to zero.
  eps_use="$eps"
  case "$label" in
    *.base) eps_use="0.0" ;;
  esac

  if ! python3 python/make_inputs.py single --kr0 "$kr0" --eps "$eps_use" \
        --cells-per-r0 "$cells" --rmax "$rmax" --n-periods "$nper" \
        --snapshots "$snaps" --outdir "$DECKDIR" >/dev/null 2>&1; then
    log "$label: FAILED to generate deck"
    printf '%s\t%s\t%s\t%s\t%s\tdeck_failed\t-\t-\t-\n' \
      "$label" "$kr0" "$eps" "$cells" "$rmax" >> "$SUMMARY"
    continue
  fi

  # make_inputs.py names decks by wavenumber and amplitude in '%g' form;
  # select the twin deck for a .base label and the perturbed deck otherwise.
  ktag=$(python3 -c 'import sys; print(f"{float(sys.argv[1]):g}")' "$kr0")
  case "$label" in
    *.base) deck="$DECKDIR/inputs.k${ktag}.base" ;;
    *)      etag=$(python3 -c 'import sys; print(f"{float(sys.argv[1]):g}")' \
                   "$eps_use")
            deck="$DECKDIR/inputs.k${ktag}.eps${etag}" ;;
  esac
  if [[ ! -f "$deck" ]]; then
    log "$label: FAILED, expected deck '$deck' was not generated"
    printf '%s\t%s\t%s\t%s\t%s\tdeck_missing\t-\t-\t-\n' \
      "$label" "$kr0" "$eps" "$cells" "$rmax" >> "$SUMMARY"
    continue
  fi

  steps=$(advance "$RUNROOT/$label" "$deck" "$label")

  # The wavelength is one axial box length and is read back from the deck so
  # that it cannot drift from the value the solver actually used.
  lambda=$(awk '/^ns.yblob/{print $3; exit}' "$RUNROOT/$label/inputs")

  rows=0
  if python3 python/dump_amplitude.py "$RUNROOT/$label" "$lambda" "$record" \
       >/dev/null 2>&1
  then
    # grep -c prints 0 and also exits non-zero when nothing matches, so the count is taken first and defaulted afterwards rather than through '|| echo 0', which would append a second line and make the arithmetic below malformed.
    rows=$(grep -cve '^[[:space:]]*#' -e '^[[:space:]]*$' "$record" 2>/dev/null || true)
    rows=${rows:-0}
  fi

  # A record with fewer than two rows cannot yield a growth rate, and the extractor exits zero when it writes a header and nothing else. Judging the case on the number of rows rather than on that exit status is what stops a run in which the solver never advanced from being recorded as a success.
  if (( rows >= 2 )); then
    tfin=$(tail -1 "$record" | awk '{print $1}')
    log "$label: complete, $steps steps, $rows snapshots, t=$tfin, record $record"

    # The figure data is exported here, while the plotfiles are still present, because this is the only moment at which they are guaranteed to exist. A failure to export is logged and does not fail the case: the growth rate is already secured by the amplitude record, and losing the interface figure of one case is not a reason to mark a completed calculation as failed.
    casebundle="$BUNDLEDIR/$label"
    mkdir -p "$casebundle"
    cp "$record" "$casebundle/amplitude.dat"
    fieldargs=()
    for f in $FIELDS; do fieldargs+=(--field "$f"); done
    if ! python3 python/dump_profiles.py "$RUNROOT/$label" "$lambda" "$casebundle" "${fieldargs[@]}" >/dev/null 2>&1; then
      log "$label: figure data export failed; amplitude record is unaffected"
    fi
    printf '%s\t%s\t%s\t%s\t%s\tok\t%s\t%s\t%s\n' \
      "$label" "$kr0" "$eps" "$cells" "$rmax" "$steps" "$tfin" "$record" \
      >> "$SUMMARY"
  else
    rm -f "$record"
    log "$label: no usable amplitude record ($rows row(s) from $steps step(s)); see $RUNROOT/$label/solver.err"
    printf '%s\t%s\t%s\t%s\t%s\tno_record\t%s\t-\t-\n' \
      "$label" "$kr0" "$eps" "$cells" "$rmax" "$steps" >> "$SUMMARY"
  fi
done < "$CASEFILE"

failed=$(awk -F'\t' 'NR>1 && $6 != "ok" { n++ } END { print n+0 }' "$SUMMARY")
if (( failed > 0 )); then
  log "sweep finished with $failed failed case(s); summary in $SUMMARY"
else
  log "sweep finished; summary in $SUMMARY"
fi
echo
if command -v column >/dev/null 2>&1; then
  column -t -s $'\t' "$SUMMARY"
else
  awk -F'\t' '{printf "%-14s %-6s %-8s %-6s %-6s %-12s %-8s %-10s\n", \
    $1, $2, $3, $4, $5, $6, $7, $8}' "$SUMMARY"
fi

# A non-zero status when any case failed, so that a scheduler records the job as failed rather than as completed. A sweep that produced no usable record is not a success, and reporting it as one is how a broken run goes unnoticed.
exit $(( failed > 0 ? 1 : 0 ))
