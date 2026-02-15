#!/usr/bin/env bash
#SBATCH --job-name=BHMimp-dmrg
#SBATCH --partition=public
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=4G
#SBATCH --requeue
#SBATCH --output=slurm-%A_%a.out

# CAIUS defaults above target full-node packing on CPU nodes.
# --mem=0 requests all memory on the allocated node.
# Override at submit time when needed, e.g.:
#   sbatch -p grant -t 4-00:00:00 --mem=64G --array=1-N hpc_campaigns/slurm/job.slurm <run_root> [threads_per_run]

set -euo pipefail # Exit on error, treat unset variables as errors, and fail on pipeline errors.

RUN_ROOT="${1:?Usage: sbatch --array=1-N hpc_campaigns/slurm/job.slurm <campaign_run_root> [threads_per_run]}" # The directory containing the jobfile and campaign data for this run.
THREADS_PER_RUN="${2:-1}" 
JOBFILE="${RUN_ROOT}/jobfile"

if [[ ! -f "$JOBFILE" ]]; then
  echo "Jobfile not found: $JOBFILE" >&2
  exit 1
fi

if [[ ! "$THREADS_PER_RUN" =~ ^[1-9][0-9]*$ ]]; then
  echo "threads_per_run must be a positive integer, got: $THREADS_PER_RUN" >&2
  exit 1
fi

REPO_ROOT="${BHM_REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
if [[ ! -f "${REPO_ROOT}/Project.toml" ]]; then
  echo "Unable to determine repository root (missing Project.toml): ${REPO_ROOT}" >&2
  echo "Set BHM_REPO_ROOT via submit_multilauncher.sh or submit with --export=ALL,BHM_REPO_ROOT=/path/to/repo" >&2
  exit 1
fi
cd "$REPO_ROOT"

LAUNCHER_SCRIPT="${REPO_ROOT}/hpc_campaigns/slurm/dynamic_multilauncher.sh"
if [[ ! -f "$LAUNCHER_SCRIPT" ]]; then
  echo "Dynamic launcher not found: $LAUNCHER_SCRIPT" >&2
  exit 1
fi

# Ensure campaign commands use this repository's Julia environment by default.
export JULIA_PROJECT="${JULIA_PROJECT:-$REPO_ROOT}"

echo "[$(date)] Starting multi-launch task"
echo "run_root=$RUN_ROOT"
echo "jobfile=$JOBFILE"
echo "array_task=${SLURM_ARRAY_TASK_ID:-unset}"
echo "partition=${SLURM_JOB_PARTITION:-unset}"
echo "time_limit_minutes=${SLURM_TIMELIMIT:-unset}"
echo "cpus_per_task=${SLURM_CPUS_PER_TASK:-unset}"
echo "mem_per_node_mb=${SLURM_MEM_PER_NODE:-unset}"
echo "threads_per_run=$THREADS_PER_RUN"
echo "julia_project=${JULIA_PROJECT}"

if [[ -n "${SLURM_CPUS_ON_NODE:-}" ]] && (( THREADS_PER_RUN > SLURM_CPUS_ON_NODE )); then
  echo "Warning: threads_per_run ($THREADS_PER_RUN) > cpus_on_node (${SLURM_CPUS_ON_NODE}); runs may serialize." >&2
fi

bash "$LAUNCHER_SCRIPT" "$JOBFILE" "$THREADS_PER_RUN"

echo "[$(date)] Finished multi-launch task"
