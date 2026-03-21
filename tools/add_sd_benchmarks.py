#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Fetch SymbolicData IntPS XML, write input/<stem>.txt, compute gb/<stem>.txt via test-f4.

Usage (from repo root):
  python3 tools/add_sd_benchmarks.py
  TEST_F4_BIN=build/test-f4 TIMEOUT=600 python3 tools/add_sd_benchmarks.py
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SD_BASE = (
    "https://raw.githubusercontent.com/symbolicdata/data/master/XMLResources/IntPS/"
)

# (output stem, IntPS filename, timeout seconds)
BENCHMARKS: list[tuple[str, str, int]] = [
    ("caprasse", "Caprasse.xml", 1800),
    ("butcher", "Butcher.xml", 600),
    ("fateman", "Fateman.xml", 600),
    ("trinks", "Trinks.xml", 300),
    ("noonburg", "Noonburg-89.xml", 600),
    ("four_body", "FourBodyProblem.xml", 1200),
    ("discriminant_4", "Discriminant_4.xml", 600),
    ("reimer_5", "Reimer_5.xml", 600),
    ("gerdt_91a", "Gerdt-91a.xml", 600),
    ("ellipsoid_2", "Ellipsoid_2.xml", 600),
]


def main() -> int:
    os.chdir(ROOT)
    py = sys.executable
    conv = ROOT / "tools" / "sd_to_pgbc_input.py"
    test_f4 = Path(os.environ.get("TEST_F4_BIN", "test/test-f4.bin"))
    if not test_f4.is_file():
        print(f"Missing {test_f4} — build with: make test/test-f4.bin", file=sys.stderr)
        return 1
    default_timeout = int(os.environ.get("TIMEOUT", "600"))

    failed: list[str] = []
    for stem, xml_name, tmo in BENCHMARKS:
        timeout = tmo if tmo else default_timeout
        url = SD_BASE + xml_name
        inp = ROOT / "input" / f"{stem}.txt"
        gb = ROOT / "gb" / f"{stem}.txt"
        print(f"=== {stem} ({xml_name}) timeout={timeout}s ===")
        try:
            subprocess.run(
                [py, str(conv), url, "-q", "-o", str(inp)],
                check=True,
                timeout=120,
            )
        except subprocess.CalledProcessError as e:
            print(f"  convert failed: {e}", file=sys.stderr)
            failed.append(stem)
            continue
        except subprocess.TimeoutExpired:
            print(f"  convert timeout", file=sys.stderr)
            failed.append(stem)
            continue
        try:
            proc = subprocess.run(
                [str(test_f4), str(inp), "1", "0", "1"],
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired:
            print(f"  test-f4 timeout after {timeout}s", file=sys.stderr)
            failed.append(stem)
            continue
        if proc.returncode != 0:
            print(proc.stderr or proc.stdout or "(no output)", file=sys.stderr)
            print(f"  test-f4 exit {proc.returncode}", file=sys.stderr)
            failed.append(stem)
            continue
        line = proc.stdout.strip()
        if not line:
            print(f"  empty GB output", file=sys.stderr)
            failed.append(stem)
            continue
        gb.write_text(line + "\n", encoding="utf-8")
        print(f"  wrote {inp.name} and {gb.name}")

    if failed:
        print(f"\nFailed: {', '.join(failed)}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
