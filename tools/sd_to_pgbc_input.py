#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Convert SymbolicData polynomial-system XML (INTPS / ModPS) to a single-line
# parallelGBC input file (same format as test/input/*.txt).
#
# Data source: https://github.com/symbolicdata/data — see tools/README-SymbolicData.md
import argparse
import re
import sys
import urllib.request
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional, Tuple


def _parse_gf_prime(text: str) -> Optional[int]:
    m = re.search(r"GF\s*\(\s*(\d+)\s*\)", text.strip(), re.I)
    return int(m.group(1)) if m else None


def _build_var_map(names: List[str]) -> Dict[str, int]:
    """Map each declared name to 1-based index x[1]..x[n] in declaration order."""
    return {name: i + 1 for i, name in enumerate(names)}


def _substitute_vars(poly: str, var_to_index: Dict[str, int]) -> str:
    """Replace algebraic variable names with x[k] tokens (parallelGBC syntax)."""
    out = poly.strip()
    # Placeholders first so a variable named 'x' does not corrupt already-emitted x[i] tokens.
    names_by_len = sorted(var_to_index.keys(), key=len, reverse=True)
    for name in names_by_len:
        idx = var_to_index[name]
        ph = f"__PGSDV{idx}__"
        pat = r"(?<![A-Za-z0-9_])" + re.escape(name) + r"(?![A-Za-z0-9_])"
        out = re.sub(pat, ph, out)
    for idx in sorted(set(var_to_index.values())):
        out = out.replace(f"__PGSDV{idx}__", f"x[{idx}]")
    return out


def _normalize_whitespace(s: str) -> str:
    return re.sub(r"\s+", "", s)


def convert_xml_bytes(data: bytes, source_hint: str = "") -> Tuple[str, Dict]:
    """
    Returns (one_line_input, meta) where meta may include 'gf_prime', 'kind', 'warnings'.
    """
    meta: Dict = {"warnings": [], "source": source_hint}
    root = ET.fromstring(data)
    tag = root.tag.split("}")[-1]  # strip namespace if present

    if tag == "INTPS":
        meta["kind"] = "INTPS"
        gf_prime = None
    elif tag == "ModPS":
        meta["kind"] = "ModPS"
        bd = root.find("basedomain")
        if bd is None or not bd.text:
            raise ValueError("ModPS: missing <basedomain>")
        gf_prime = _parse_gf_prime(bd.text)
        if gf_prime is None:
            raise ValueError(f"ModPS: could not parse GF prime from {bd.text!r}")
        meta["gf_prime"] = gf_prime
        if gf_prime != 32003:
            meta["warnings"].append(
                f"GF({gf_prime}) != GF(32003): test-f4 uses 32003; coefficients may not match after bringIn."
            )
    else:
        raise ValueError(f"Unsupported root element <{tag}> (expected INTPS or ModPS)")

    vel = root.find("vars")
    if vel is None or not vel.text:
        raise ValueError("Missing <vars>")
    raw_vars = [v.strip() for v in vel.text.split(",") if v.strip()]
    if not raw_vars:
        raise ValueError("Empty <vars>")

    bel = root.find("basis")
    if bel is None:
        raise ValueError("Missing <basis>")
    polys = []
    for p in bel.findall("poly"):
        if p.text is None:
            continue
        t = p.text.strip()
        if t:
            polys.append(t)
    if not polys:
        raise ValueError("No <poly> entries under <basis>")

    vmap = _build_var_map(raw_vars)
    mapped = [_normalize_whitespace(_substitute_vars(poly, vmap)) for poly in polys]
    line = ", ".join(mapped)
    meta["n_vars"] = len(raw_vars)
    meta["n_polys"] = len(mapped)
    return line, meta


def _read_input(path_or_url: str) -> Tuple[bytes, str]:
    if re.match(r"^https?://", path_or_url, re.I):
        req = urllib.request.Request(
            path_or_url,
            headers={"User-Agent": "parallelGBC-sd_to_pgbc_input/1.0"},
        )
        with urllib.request.urlopen(req, timeout=60) as resp:
            return resp.read(), path_or_url
    with open(path_or_url, "rb") as f:
        return f.read(), path_or_url


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Convert SymbolicData INTPS/ModPS XML to one-line parallelGBC input."
    )
    ap.add_argument(
        "input",
        help="Path to .xml or http(s) URL (e.g. raw.githubusercontent.com/.../Cyclic_4.xml)",
    )
    ap.add_argument(
        "-o",
        "--output",
        help="Write here instead of stdout",
    )
    ap.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress stderr metadata / warnings",
    )
    args = ap.parse_args()

    try:
        raw, hint = _read_input(args.input)
        line, meta = convert_xml_bytes(raw, source_hint=hint)
    except (ET.ParseError, ValueError, OSError, urllib.error.URLError) as e:
        print(f"sd_to_pgbc_input: {e}", file=sys.stderr)
        return 1

    if not args.quiet:
        print(
            f"# kind={meta['kind']} n_vars={meta['n_vars']} n_polys={meta['n_polys']}",
            file=sys.stderr,
        )
        if meta.get("gf_prime") is not None:
            print(f"# gf_prime={meta['gf_prime']}", file=sys.stderr)
        for w in meta.get("warnings", []):
            print(f"# WARNING: {w}", file=sys.stderr)

    if args.output:
        with open(args.output, "w", encoding="utf-8") as out:
            out.write(line + "\n")
    else:
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
