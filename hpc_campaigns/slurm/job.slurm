#!/usr/bin/env bash
#SBATCH --job-name=bhm-dmrg
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=slurm-%A_%a.out

set -euo pipefail

RUN_ROOT="${1:?Usage: sbatch hpc_campaigns/slurm/job.slurm <campaign_run_root>}"
INDEX_FILE="${RUN_ROOT}/index.csv"
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is not set}"

if [[ ! -f "$INDEX_FILE" ]]; then
  echo "Index file not found: $INDEX_FILE" >&2
  exit 1
fi

LINE="$(sed -n "$((TASK_ID + 2))p" "$INDEX_FILE")"
if [[ -z "$LINE" ]]; then
  echo "No run entry for array task $TASK_ID" >&2
  exit 1
fi

IFS=',' read -r RUN_ID RUN_DIR PARAMS_PATH APP_SCRIPT STATUS <<< "$LINE"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

echo "[$(date)] Starting $RUN_ID"
echo "run_dir=$RUN_DIR"
echo "params=$PARAMS_PATH"
echo "app=$APP_SCRIPT"

julia "$APP_SCRIPT" "$PARAMS_PATH"

echo "[$(date)] Finished $RUN_ID"
