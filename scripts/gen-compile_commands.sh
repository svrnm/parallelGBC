#!/usr/bin/env bash
# Generate compile_commands.json at repo root for clangd / IDE.
# Requires: https://github.com/rizsotto/Bear (brew install bear, apt install bear, etc.)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
if ! command -v bear >/dev/null 2>&1; then
	echo "Install Bear (e.g. brew install bear) and re-run." >&2
	exit 1
fi
bear -- make clean all
echo "Wrote $ROOT/compile_commands.json"
