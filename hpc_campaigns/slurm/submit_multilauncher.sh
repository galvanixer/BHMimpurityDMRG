#!/usr/bin/env bash
set -euo pipefail

# Submit a campaign "jobfile" through Slurm array tasks.
# Each array task runs job.slurm, which calls:
#   dynamic_multilauncher.sh <run_root>/jobfile <threads_per_run>
#
# Args:
#   1) run_root        : campaign directory containing "jobfile" (required)
#   2) num_array_tasks : Slurm array task count; omit or set 0 for auto
#                        auto = ceil(jobfile_lines / 24) using 24-core standard nodes
#                        On CAIUS full-node partitions this is effectively the node count.
#   3) threads_per_run : CPU threads used by each launched command inside a task
#   4) cpus_per_task   : Slurm CPU request per array task (default: 24)
#   5) partition       : optional submit override (e.g., public, grant, publicgpu, grantgpu)
#   6) account_name    : optional account selector for grant/grantgpu
#   7) job_script      : optional alternate Slurm worker script path
# Order:
#   <run_root> [num_array_tasks] [threads_per_run] [cpus_per_task] [partition] [account_name] [job_script]
#                        - profile key (e.g., francesco -> CAIUS_GRANT_ACCOUNT_FRANCESCO)
#                        - direct account via acct:<account_id> (e.g., acct:g2025a457b)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" # Absolute directory of this script, used for resolving relative paths reliably.
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)" # Repository root, assumed to contain Project.toml. Can be overridden by BHM_REPO_ROOT env var.
LOCAL_CREDENTIALS_FILE="${SCRIPT_DIR}/credentials.local.sh" # Optional local credentials file, kept out of git, for partition-bound accounts.
STANDARD_NODE_CORES=24 # Standard CPU cores per node on CAIUS; used for auto array task calculation. Adjust if using a different cluster with different node sizes.

RUN_ROOT="${1:?Usage: bash hpc_campaigns/slurm/submit_multilauncher.sh <campaign_run_root> [num_array_tasks] [threads_per_run] [cpus_per_task] [partition] [account_name] [job_script]}"
NUM_ARRAY_TASKS_RAW="${2:-}"
THREADS_PER_RUN="${3:-1}"
CPUS_PER_TASK_REQUEST="${4:-$STANDARD_NODE_CORES}"
PARTITION="${5:-}"
ACCOUNT_NAME="${6:-}"
JOB_SCRIPT_RAW="${7:-hpc_campaigns/slurm/job.slurm}"
JOBFILE="${RUN_ROOT}/jobfile"

# Optional local credentials (kept out of git) for partition-bound accounts.
if [[ -f "$LOCAL_CREDENTIALS_FILE" ]]; then
  # shellcheck disable=SC1090
  source "$LOCAL_CREDENTIALS_FILE"
fi

# Resolve account with priority:
# 1) selector-specific prefixed var (e.g., CAIUS_GRANT_ACCOUNT_FRANCESCO)
# 2) selector-specific legacy var   (e.g., GRANT_ACCOUNT_FRANCESCO)
# 3) user-specific prefixed var     (e.g., CAIUS_GRANT_ACCOUNT_TANUL)
# 4) user-specific legacy var       (e.g., GRANT_ACCOUNT_TANUL)
# 5) global prefixed var            (e.g., CAIUS_GRANT_ACCOUNT)
# 6) global legacy var              (e.g., GRANT_ACCOUNT)
resolve_partition_account() {
  local partition_kind="$1"
  local selector="${2:-}"
  local selector_key
  local user_key
  local v

  if [[ -n "$selector" ]]; then
    selector_key="$(printf '%s' "$selector" | tr '[:lower:]' '[:upper:]' | tr -c 'A-Z0-9' '_')"
    v="CAIUS_${partition_kind}_ACCOUNT_${selector_key}"
    if [[ -n "${!v:-}" ]]; then
      printf '%s' "${!v}"
      return
    fi

    v="${partition_kind}_ACCOUNT_${selector_key}"
    if [[ -n "${!v:-}" ]]; then
      printf '%s' "${!v}"
      return
    fi
  fi

  user_key="$(printf '%s' "${USER:-}" | tr '[:lower:]' '[:upper:]' | tr -c 'A-Z0-9' '_')"
  if [[ -n "$user_key" ]]; then
    v="CAIUS_${partition_kind}_ACCOUNT_${user_key}"
    if [[ -n "${!v:-}" ]]; then
      printf '%s' "${!v}"
      return
    fi

    v="${partition_kind}_ACCOUNT_${user_key}"
    if [[ -n "${!v:-}" ]]; then
      printf '%s' "${!v}"
      return
    fi
  fi

  v="CAIUS_${partition_kind}_ACCOUNT"
  if [[ -n "${!v:-}" ]]; then
    printf '%s' "${!v}"
    return
  fi

  v="${partition_kind}_ACCOUNT"
  if [[ -n "${!v:-}" ]]; then
    printf '%s' "${!v}"
    return
  fi
}

# Validate campaign layout and basic numeric inputs early.
if [[ ! -d "$RUN_ROOT" ]]; then
  echo "Run root directory not found: $RUN_ROOT" >&2
  exit 1
fi

if [[ ! -f "$JOBFILE" ]]; then
  echo "Missing jobfile: $JOBFILE" >&2
  exit 1
fi

RUN_ROOT="$(cd "$RUN_ROOT" && pwd)"
JOBFILE="${RUN_ROOT}/jobfile"

if [[ "$JOB_SCRIPT_RAW" = /* ]]; then
  JOB_SCRIPT="$JOB_SCRIPT_RAW"
elif [[ -f "$JOB_SCRIPT_RAW" ]]; then
  JOB_SCRIPT="$(cd "$(dirname "$JOB_SCRIPT_RAW")" && pwd)/$(basename "$JOB_SCRIPT_RAW")"
elif [[ -f "${REPO_ROOT}/${JOB_SCRIPT_RAW}" ]]; then
  JOB_SCRIPT="${REPO_ROOT}/${JOB_SCRIPT_RAW}"
else
  echo "Unable to locate job script: $JOB_SCRIPT_RAW" >&2
  echo "Tried current directory and repository root (${REPO_ROOT})." >&2
  exit 1
fi

if [[ ! "$THREADS_PER_RUN" =~ ^[1-9][0-9]*$ ]]; then
  echo "threads_per_run must be a positive integer, got: $THREADS_PER_RUN" >&2
  exit 1
fi

if [[ ! "$CPUS_PER_TASK_REQUEST" =~ ^[1-9][0-9]*$ ]]; then
  echo "cpus_per_task must be a positive integer, got: $CPUS_PER_TASK_REQUEST" >&2
  exit 1
fi

# Number of commands in jobfile bounds useful array size.
N_CMDS="$(awk 'END { print NR }' "$JOBFILE")"
if (( N_CMDS <= 0 )); then
  echo "No commands found in $JOBFILE" >&2
  exit 1
fi

if [[ -n "$NUM_ARRAY_TASKS_RAW" ]] && [[ ! "$NUM_ARRAY_TASKS_RAW" =~ ^[0-9]+$ ]]; then
  echo "num_array_tasks must be a non-negative integer, got: $NUM_ARRAY_TASKS_RAW" >&2
  exit 1
fi

# num_array_tasks policy:
# - omitted or 0: auto = ceil(N_CMDS / 24), using 24-core standard node assumption.
# - > N_CMDS: clamp to N_CMDS (never schedule more tasks than commands).
if [[ -z "$NUM_ARRAY_TASKS_RAW" || "$NUM_ARRAY_TASKS_RAW" == "0" ]]; then
  NUM_ARRAY_TASKS="$(( (N_CMDS + STANDARD_NODE_CORES - 1) / STANDARD_NODE_CORES ))"
  ARRAY_TASKS_SOURCE="auto ceil(commands/${STANDARD_NODE_CORES})"
else
  NUM_ARRAY_TASKS="$NUM_ARRAY_TASKS_RAW"
  ARRAY_TASKS_SOURCE="user"
fi

if (( NUM_ARRAY_TASKS > N_CMDS )); then
  NUM_ARRAY_TASKS=$N_CMDS
fi

# Submit a 1-based contiguous array range; job.slurm uses SLURM_ARRAY_TASK_ID.
# On CAIUS, each array task usually maps to one full node, so num_array_tasks ~= nodes.
ARRAY_SPEC="1-${NUM_ARRAY_TASKS}"
SBATCH_ARGS=(--array="$ARRAY_SPEC" --cpus-per-task="$CPUS_PER_TASK_REQUEST" --export="ALL,BHM_REPO_ROOT=${REPO_ROOT}")
ACCOUNT_TO_USE=""
ACCOUNT_ENV_HINT=""
PARTITION_KIND=""
TIME_LIMIT_OVERRIDE=""

if [[ -n "$PARTITION" ]]; then
  if [[ ! "$PARTITION" =~ ^[A-Za-z0-9_-]+$ ]]; then
    echo "Invalid partition: $PARTITION" >&2
    exit 1
  fi
  SBATCH_ARGS+=(--partition="$PARTITION")

  case "$PARTITION" in
    grant)
      PARTITION_KIND="GRANT"
      if [[ "$ACCOUNT_NAME" == acct:* ]]; then
        ACCOUNT_TO_USE="${ACCOUNT_NAME#acct:}"
      else
        ACCOUNT_TO_USE="$(resolve_partition_account "GRANT" "$ACCOUNT_NAME")"
      fi
      ACCOUNT_ENV_HINT="CAIUS_GRANT_ACCOUNT"
      ;;
    grantgpu)
      PARTITION_KIND="GRANTGPU"
      if [[ "$ACCOUNT_NAME" == acct:* ]]; then
        ACCOUNT_TO_USE="${ACCOUNT_NAME#acct:}"
      else
        ACCOUNT_TO_USE="$(resolve_partition_account "GRANTGPU" "$ACCOUNT_NAME")"
      fi
      ACCOUNT_ENV_HINT="CAIUS_GRANTGPU_ACCOUNT"
      ;;
    *)
      ACCOUNT_TO_USE=""
      ;;
  esac

  if [[ "$PARTITION" == "grant" || "$PARTITION" == "grantgpu" ]]; then
    if [[ -z "$ACCOUNT_TO_USE" ]]; then
      echo "Missing account credential for partition '$PARTITION'." >&2
      if [[ -n "$ACCOUNT_NAME" ]]; then
        echo "No account found for selector '$ACCOUNT_NAME'." >&2
        echo "Define ${ACCOUNT_ENV_HINT}_<NAME> (or ${PARTITION_KIND}_ACCOUNT_<NAME>) in ${LOCAL_CREDENTIALS_FILE}," >&2
        echo "or pass direct account as acct:<account_id>." >&2
      else
        echo "Set env var ${ACCOUNT_ENV_HINT} or define it in ${LOCAL_CREDENTIALS_FILE}." >&2
      fi
      echo "Also supports user-specific vars like CAIUS_GRANT_ACCOUNT_<USER_KEY> / CAIUS_GRANTGPU_ACCOUNT_<USER_KEY>." >&2
      exit 1
    fi
    SBATCH_ARGS+=(--account="$ACCOUNT_TO_USE")
  fi
fi

# Partition-based walltime defaults. Submit-time --time overrides job.slurm header.
case "${PARTITION:-public}" in
  public|publicgpu)
    TIME_LIMIT_OVERRIDE="1-00:00:00"
    ;;
  grant|grantgpu)
    TIME_LIMIT_OVERRIDE="4-00:00:00"
    ;;
  *)
    TIME_LIMIT_OVERRIDE=""
    ;;
esac

if [[ -n "$TIME_LIMIT_OVERRIDE" ]]; then
  SBATCH_ARGS+=(--time="$TIME_LIMIT_OVERRIDE")
fi

echo "Submitting multilaunch campaign"
echo "commands:         $N_CMDS"
echo "num_array_tasks:  $NUM_ARRAY_TASKS"
echo "array_tasks_src:  $ARRAY_TASKS_SOURCE"
echo "threads_per_run:  $THREADS_PER_RUN"
echo "array spec:       $ARRAY_SPEC"
echo "job script:       $JOB_SCRIPT"
echo "run root:         $RUN_ROOT"
echo "repo root:        $REPO_ROOT"
echo "partition:        ${PARTITION:-<job.slurm default>}"
echo "time_limit:       ${TIME_LIMIT_OVERRIDE:-<job.slurm default>}"
echo "cpus_per_task:    $CPUS_PER_TASK_REQUEST"
echo "account_name:     ${ACCOUNT_NAME:-<auto>}"
if [[ -n "$ACCOUNT_TO_USE" ]]; then
  echo "account:          $ACCOUNT_TO_USE"
fi
if (( NUM_ARRAY_TASKS > 1 )); then
  echo "NOTE: this cluster may allocate one full node per array task; num_array_tasks=$NUM_ARRAY_TASKS can reserve multiple full nodes"
fi

# Forward run_root and threads_per_run as positional args to job.slurm.
sbatch "${SBATCH_ARGS[@]}" "$JOB_SCRIPT" "$RUN_ROOT" "$THREADS_PER_RUN"
