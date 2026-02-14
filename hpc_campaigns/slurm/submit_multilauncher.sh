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
#   5) partition       : optional submit override (e.g., public, grant, publicgpu, grantgpu)
#   6) account_name    : optional account selector for grant/grantgpu
#                        - profile key (e.g., francesco -> CAIUS_GRANT_ACCOUNT_FRANCESCO)
#                        - direct account via acct:<account_id> (e.g., acct:g2025a457b)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOCAL_CREDENTIALS_FILE="${SCRIPT_DIR}/credentials.local.sh"

RUN_ROOT="${1:?Usage: bash hpc_campaigns/slurm/submit_multilauncher.sh <campaign_run_root> [n_subjobs] [threads_per_run] [job_script] [partition] [account_name]}"
N_SUBJOBS="${2:-0}"
THREADS_PER_RUN="${3:-1}"
JOB_SCRIPT="${4:-hpc_campaigns/slurm/job.slurm}"
PARTITION="${5:-}"
ACCOUNT_NAME="${6:-}"
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
SBATCH_ARGS=(--array="$ARRAY_SPEC")
ACCOUNT_TO_USE=""
ACCOUNT_ENV_HINT=""
PARTITION_KIND=""

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

echo "Submitting multilaunch campaign"
echo "commands:         $N_CMDS"
echo "subjobs (array):  $N_SUBJOBS"
echo "threads_per_run:  $THREADS_PER_RUN"
echo "array spec:       $ARRAY_SPEC"
echo "job script:       $JOB_SCRIPT"
echo "run root:         $RUN_ROOT"
echo "partition:        ${PARTITION:-<job.slurm default>}"
echo "account_name:     ${ACCOUNT_NAME:-<auto>}"
if [[ -n "$ACCOUNT_TO_USE" ]]; then
  echo "account:          $ACCOUNT_TO_USE"
fi

# Forward run_root and threads_per_run as positional args to job.slurm.
sbatch "${SBATCH_ARGS[@]}" "$JOB_SCRIPT" "$RUN_ROOT" "$THREADS_PER_RUN"
