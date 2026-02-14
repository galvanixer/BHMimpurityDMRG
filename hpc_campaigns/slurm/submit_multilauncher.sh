#!/usr/bin/env bash
set -euo pipefail

# Submit a campaign "jobfile" through Slurm array tasks.
# Each array task runs job.slurm, which calls:
#   hpc_multilauncher <run_root>/jobfile <threads_per_run>
#
# Args:
#   1) run_root        : campaign directory containing "jobfile" (required)
#   2) n_subjobs       : Slurm array task count; 0 means "auto = one task per command"
#   3) threads_per_run : CPU threads used by each launched command inside a task
#   4) job_script      : Slurm worker script path (defaults to hpc_campaigns/slurm/job.slurm)

RUN_ROOT="${1:?Usage: bash hpc_campaigns/slurm/submit_multilauncher.sh <campaign_run_root> [n_subjobs] [threads_per_run] [job_script]}"
N_SUBJOBS="${2:-0}"
THREADS_PER_RUN="${3:-1}"
JOB_SCRIPT="${4:-hpc_campaigns/slurm/job.slurm}"
JOBFILE="${RUN_ROOT}/jobfile"

# Validate campaign layout and basic numeric inputs early.
if [[ ! -d "$RUN_ROOT" ]]; then
  echo "Run root directory not found: $RUN_ROOT" >&2
  exit 1
fi

if [[ ! -f "$JOBFILE" ]]; then
  echo "Missing jobfile: $JOBFILE" >&2
  exit 1
fi

if [[ ! "$THREADS_PER_RUN" =~ ^[1-9][0-9]*$ ]]; then
  echo "threads_per_run must be a positive integer, got: $THREADS_PER_RUN" >&2
  exit 1
fi

# Number of commands in jobfile bounds useful array size.
N_CMDS=$(wc -l < "$JOBFILE")
if (( N_CMDS <= 0 )); then
  echo "No commands found in $JOBFILE" >&2
  exit 1
fi

if [[ ! "$N_SUBJOBS" =~ ^[0-9]+$ ]]; then
  echo "n_subjobs must be a non-negative integer, got: $N_SUBJOBS" >&2
  exit 1
fi

# n_subjobs policy:
# - 0: auto-expand to one array task per command.
# - > N_CMDS: clamp to N_CMDS (never schedule more tasks than commands).
if (( N_SUBJOBS == 0 )); then
  N_SUBJOBS=$N_CMDS
fi

if (( N_SUBJOBS > N_CMDS )); then
  N_SUBJOBS=$N_CMDS
fi

# Submit a 1-based contiguous array range; job.slurm uses SLURM_ARRAY_TASK_ID.
ARRAY_SPEC="1-${N_SUBJOBS}"
echo "Submitting multilaunch campaign"
echo "commands:         $N_CMDS"
echo "subjobs (array):  $N_SUBJOBS"
echo "threads_per_run:  $THREADS_PER_RUN"
echo "array spec:       $ARRAY_SPEC"
echo "job script:       $JOB_SCRIPT"
echo "run root:         $RUN_ROOT"

# Forward run_root and threads_per_run as positional args to job.slurm.
sbatch --array="$ARRAY_SPEC" "$JOB_SCRIPT" "$RUN_ROOT" "$THREADS_PER_RUN"
