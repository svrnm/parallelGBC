#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Compare committed input/*.txt to fresh SymbolicData XML conversion (multiset of generators)."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Repo root (parent of tools/)
ROOT = Path(__file__).resolve().parent.parent
SD_BASE = (
    "https://raw.githubusercontent.com/symbolicdata/data/master/XMLResources/IntPS/"
)

# (input stem without .txt, SymbolicData XML filename under IntPS/)
# Only INTPS: integer coefficients; test-f4 reduces mod 32003.
# Pairs are limited to what exists on github.com/symbolicdata/data IntPS and matches
# the committed generator *multiset* (ordering of <poly> tags may differ).
VERIFIED_PAIRS: list[tuple[str, str]] = [
    ("cyclic4", "Cyclic_4.xml"),
    ("cyclic5", "Cyclic_5.xml"),
    ("cyclic6", "Cyclic_6.xml"),
    ("cyclic7", "Cyclic_7.xml"),
    ("cyclic9", "Cyclic_9.xml"),
    # IntPS has no Cyclic_3 / Cyclic_10; cyclic8 in-tree differs from Cyclic_8.xml generators.
    # Katsura_7.xml in SD is 8 polynomials in 8 variables; katsura7.txt here is the 7-var form.
    ("katsura8", "Katsura_8.xml"),
]


def _normalize_input_line(content: str) -> list[str]:
    """Split one-line input into generator strings; strip whitespace; multiset compare."""
    t = "".join(content.split())
    parts = [p.strip() for p in t.replace(";", ",").split(",")]
    return sorted(p for p in parts if p)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Verify sd_to_pgbc_input output matches committed input (generator multiset)."
    )
    ap.add_argument(
        "--root",
        type=Path,
        default=ROOT,
        help="parallelGBC root (default: auto)",
    )
    args = ap.parse_args()
    root = args.root
    sys.path.insert(0, str(root / "tools"))
    from sd_to_pgbc_input import convert_xml_bytes  # noqa: E402
    import urllib.request

    failed = 0
    for stem, xml_name in VERIFIED_PAIRS:
        inp = root / "input" / f"{stem}.txt"
        if not inp.is_file():
            print(f"SKIP {stem}: missing {inp}")
            continue
        url = SD_BASE + xml_name
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "parallelGBC-verify_sd/1.0"}
            )
            with urllib.request.urlopen(req, timeout=120) as r:
                raw = r.read()
            line, meta = convert_xml_bytes(raw, source_hint=url)
        except Exception as e:
            print(f"FAIL {stem}: fetch/convert {xml_name}: {e}")
            failed += 1
            continue
        committed = _normalize_input_line(inp.read_text(encoding="utf-8"))
        fresh = _normalize_input_line(line)
        if committed == fresh:
            print(f"OK   {stem} ({len(committed)} generators)")
        else:
            print(f"FAIL {stem}: multiset mismatch vs {xml_name}")
            failed += 1
            if len(committed) != len(fresh):
                print(f"      counts: committed={len(committed)} sd={len(fresh)}")
            else:
                for i, (a, b) in enumerate(zip(committed, fresh)):
                    if a != b:
                        print(f"      first diff at sorted index {i}:")
                        print(f"      A: {a[:120]}...")
                        print(f"      B: {b[:120]}...")
                        break

    if failed:
        print(f"\n{failed} pair(s) failed.", file=sys.stderr)
        return 1
    print("\nAll available pairs matched.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
