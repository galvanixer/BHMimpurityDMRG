#!/usr/bin/env bash
#SBATCH --job-name=bhm-dmrg
#SBATCH --partition=public
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --requeue
#SBATCH --output=slurm-%A_%a.out

# CAIUS defaults above target fast scheduling on CPU nodes.
# Override at submit time when needed, e.g.:
#   sbatch -p grant -t 4-00:00:00 --mem=64G --array=1-N hpc_campaigns/slurm/job.slurm <run_root> [threads_per_run]

set -euo pipefail

RUN_ROOT="${1:?Usage: sbatch --array=1-N hpc_campaigns/slurm/job.slurm <campaign_run_root> [threads_per_run]}"
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

if ! command -v hpc_multilauncher >/dev/null 2>&1; then
  echo "hpc_multilauncher not found in PATH" >&2
  exit 1
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

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

if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]] && (( THREADS_PER_RUN > SLURM_CPUS_PER_TASK )); then
  echo "Warning: threads_per_run ($THREADS_PER_RUN) > cpus_per_task (${SLURM_CPUS_PER_TASK}); runs may serialize." >&2
fi

hpc_multilauncher "$JOBFILE" "$THREADS_PER_RUN"

echo "[$(date)] Finished multi-launch task"
