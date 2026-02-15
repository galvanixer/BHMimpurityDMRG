#!/usr/bin/env bash
set -euo pipefail

COMMAND_FILE="${1:?Usage: dynamic_multilauncher.sh <command_file> [threads_per_run]}"
THREADS_PER_RUN="${2:-1}"

if [[ ! -f "$COMMAND_FILE" ]]; then
  echo "Command file not found: $COMMAND_FILE" >&2
  exit 1
fi

if [[ ! "$THREADS_PER_RUN" =~ ^[1-9][0-9]*$ ]]; then
  echo "threads_per_run must be a positive integer, got: $THREADS_PER_RUN" >&2
  exit 1
fi

TOTAL_CMDS="$(awk 'END { print NR }' "$COMMAND_FILE")"
if (( TOTAL_CMDS <= 0 )); then
  echo "No commands found in $COMMAND_FILE"
  exit 0
fi

RUN_ROOT="$(cd "$(dirname "$COMMAND_FILE")" && pwd)"
JOB_GROUP_ID="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
STATE_DIR="${RUN_ROOT}/.dynamic_launcher_state/${JOB_GROUP_ID}"
LOCK_DIR="${STATE_DIR}/lock"
NEXT_FILE="${STATE_DIR}/next_line"
FAIL_COUNT_FILE="${STATE_DIR}/fail_count"
FAIL_LOG_FILE="${STATE_DIR}/failed_commands.tsv"

mkdir -p "$STATE_DIR"

acquire_lock() {
  while ! mkdir "$LOCK_DIR" 2>/dev/null; do
    sleep 0.05
  done
}

release_lock() {
  rmdir "$LOCK_DIR" 2>/dev/null || true
}

# Initialize shared state once for the whole array job.
acquire_lock
if [[ ! -f "$NEXT_FILE" ]]; then
  echo 1 > "$NEXT_FILE"
fi
if [[ ! -f "$FAIL_COUNT_FILE" ]]; then
  echo 0 > "$FAIL_COUNT_FILE"
fi
release_lock

record_failure() {
  local line_no="$1"
  local rc="$2"
  local cmd="$3"
  acquire_lock
  local n_fail
  n_fail="$(cat "$FAIL_COUNT_FILE")"
  echo $((n_fail + 1)) > "$FAIL_COUNT_FILE"
  printf "%s\t%s\t%s\n" "$line_no" "$rc" "$cmd" >> "$FAIL_LOG_FILE"
  release_lock
}

get_next_line_number() {
  local next
  acquire_lock
  next="$(cat "$NEXT_FILE")"
  if (( next > TOTAL_CMDS )); then
    release_lock
    return 1
  fi
  echo $((next + 1)) > "$NEXT_FILE"
  release_lock
  printf "%s\n" "$next"
}

CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-}"
CPUS_ON_NODE="${SLURM_CPUS_ON_NODE:-}"

if [[ -n "$CPUS_PER_TASK" && -n "$CPUS_ON_NODE" ]]; then
  if (( CPUS_PER_TASK < CPUS_ON_NODE )); then
    CPU_BUDGET="$CPUS_PER_TASK"
  else
    CPU_BUDGET="$CPUS_ON_NODE"
  fi
elif [[ -n "$CPUS_PER_TASK" ]]; then
  CPU_BUDGET="$CPUS_PER_TASK"
elif [[ -n "$CPUS_ON_NODE" ]]; then
  CPU_BUDGET="$CPUS_ON_NODE"
else
  CPU_BUDGET=1
fi

PARALLEL_SLOTS=$((CPU_BUDGET / THREADS_PER_RUN))
if (( PARALLEL_SLOTS < 1 )); then
  PARALLEL_SLOTS=1
fi

echo "dynamic_multilauncher: commands=${TOTAL_CMDS} cpu_budget=${CPU_BUDGET} cpus_per_task=${CPUS_PER_TASK:-unset} cpus_on_node=${CPUS_ON_NODE:-unset} threads_per_run=${THREADS_PER_RUN} slots=${PARALLEL_SLOTS}"
echo "dynamic_multilauncher: state_dir=${STATE_DIR}"

worker_loop() {
  local worker_id="$1"
  local line_no cmd rc
  while true; do
    if ! line_no="$(get_next_line_number)"; then
      break
    fi
    cmd="$(sed -n "${line_no}p" "$COMMAND_FILE")"
    if [[ -z "${cmd//[[:space:]]/}" ]]; then
      continue
    fi

    echo "[$(date)] [worker=${worker_id} line=${line_no}] START"
    if bash -lc "$cmd"; then
      echo "[$(date)] [worker=${worker_id} line=${line_no}] OK"
    else
      rc=$?
      echo "[$(date)] [worker=${worker_id} line=${line_no}] FAIL rc=${rc}" >&2
      record_failure "$line_no" "$rc" "$cmd"
    fi
  done
}

declare -a pids=()
for worker_id in $(seq 1 "$PARALLEL_SLOTS"); do
  worker_loop "$worker_id" &
  pids+=("$!")
done

for pid in "${pids[@]}"; do
  wait "$pid"
done

N_FAIL="$(cat "$FAIL_COUNT_FILE")"
if (( N_FAIL > 0 )); then
  echo "dynamic_multilauncher: ${N_FAIL} command(s) failed. See ${FAIL_LOG_FILE}" >&2
  exit 1
fi

echo "dynamic_multilauncher: all commands completed successfully"
