#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${1:-${BHM_REPO_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}}"
JULIA_DEPOT="${2:-${BHM_JULIA_DEPOT:-${JULIA_DEPOT_PATH:-$HOME/.julia_bhmimpuritydmrg}}}"

if [[ ! -f "${REPO_ROOT}/Project.toml" ]]; then
  echo "Project.toml not found under repo root: ${REPO_ROOT}" >&2
  exit 1
fi

export JULIA_PROJECT="$REPO_ROOT"
export JULIA_DEPOT_PATH="$JULIA_DEPOT"
export JULIA_CPU_TARGET="${JULIA_CPU_TARGET:-generic}"
export JULIA_PKG_PRECOMPILE_AUTO=0

DEPOT_PRIMARY="${JULIA_DEPOT_PATH%%:*}"
if [[ -n "$DEPOT_PRIMARY" ]]; then
  mkdir -p "$DEPOT_PRIMARY"
fi

echo "Precompiling Julia project..."
echo "repo_root=${REPO_ROOT}"
echo "julia_depot_path=${JULIA_DEPOT_PATH}"
echo "julia_cpu_target=${JULIA_CPU_TARGET}"

julia --project="$REPO_ROOT" -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

echo "Done."
