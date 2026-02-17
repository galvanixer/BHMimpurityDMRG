#!/usr/bin/env bash
# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  echo "Please source this file: source $(realpath "$0")" >&2
  exit 1
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
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
