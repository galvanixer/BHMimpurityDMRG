#!/usr/bin/env bash
# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

_bhm_script_path=""
_bhm_sourced=0

if [[ -n "${BASH_VERSION-}" ]]; then
  _bhm_script_path="${BASH_SOURCE[0]}"
  [[ "${BASH_SOURCE[0]}" != "${0}" ]] && _bhm_sourced=1
elif [[ -n "${ZSH_VERSION-}" ]]; then
  _bhm_script_path="${(%):-%N}"
  case "${ZSH_EVAL_CONTEXT-}" in
    *:file) _bhm_sourced=1 ;;
    *) _bhm_sourced=0 ;;
  esac
else
  _bhm_script_path="$0"
  _bhm_sourced=0
fi

if [[ $_bhm_sourced -eq 0 ]]; then
  echo "Please source this file: source ${_bhm_script_path}" >&2
  exit 1
fi

REPO_ROOT="$(cd "$(dirname "${_bhm_script_path}")" && pwd)"
BIN_DIR="$REPO_ROOT/bin"

if [[ ! -d "$BIN_DIR" ]]; then
  echo "bin directory not found: $BIN_DIR" >&2
  return 1
fi

case ":$PATH:" in
  *":$BIN_DIR:"*) ;;
  *) export PATH="$BIN_DIR:$PATH" ;;
esac

export BHMIMPURITYDMRG_ROOT="$REPO_ROOT"

unset _bhm_script_path
unset _bhm_sourced
