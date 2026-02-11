#!/usr/bin/env bash
set -euo pipefail

RUN_ROOT="${1:?Usage: bash hpc_campaigns/slurm/submit_array.sh <campaign_run_root> [job_script]}"
JOB_SCRIPT="${2:-hpc_campaigns/slurm/job.slurm}"

if [[ ! -d "$RUN_ROOT" ]]; then
  echo "Run root directory not found: $RUN_ROOT" >&2
  exit 1
fi

INDEX_FILE="${RUN_ROOT}/index.csv"
if [[ ! -f "$INDEX_FILE" ]]; then
  echo "Missing index file: $INDEX_FILE" >&2
  exit 1
fi

N_RUNS=$(( $(wc -l < "$INDEX_FILE") - 1 ))
if (( N_RUNS <= 0 )); then
  echo "No runs found in $INDEX_FILE" >&2
  exit 1
fi

ARRAY_SPEC="0-$((N_RUNS - 1))"
echo "Submitting $N_RUNS runs with array spec $ARRAY_SPEC"
echo "Job script: $JOB_SCRIPT"
echo "Run root:   $RUN_ROOT"

sbatch --array="$ARRAY_SPEC" "$JOB_SCRIPT" "$RUN_ROOT"
