#!/usr/bin/env bash
# Symlink build/compile_commands.json to the repo root for clangd / IDEs.
# Run after: cmake -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON && cmake --build build
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
if [ ! -f build/compile_commands.json ]; then
	echo "Missing build/compile_commands.json. Run:" >&2
	echo "  cmake -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Release" >&2
	echo "  cmake --build build -j\$(nproc)" >&2
	exit 1
fi
ln -sf build/compile_commands.json "$ROOT/compile_commands.json"
echo "Linked $ROOT/compile_commands.json -> build/compile_commands.json"
