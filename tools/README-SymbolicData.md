# SymbolicData → parallelGBC input

[`sd_to_pgbc_input.py`](sd_to_pgbc_input.py) turns **SymbolicData** XML (`INTPS` over **Z**, or `ModPS` over **GF(p)**) into a **single-line** file in the same format as [`test/input/`](../test/input/): variables become `x[1]…x[n]` in **declaration order** (important for degrevlex), polynomials separated by `", "`.

## License

The **SymbolicData** project states that tools and data are available under the **GNU GPL**; see [their documentation](https://symbolicdata.readthedocs.io/en/latest/). If you redistribute converted snippets or upstream XML, keep license and attribution requirements in mind. This helper script is **GPL-3.0-or-later**, consistent with parallelGBC.

## Usage

```bash
# Local clone or download
./tools/sd_to_pgbc_input.py path/to/Cyclic_4.xml -o /tmp/cyclic4_pgbc.txt

# Directly from GitHub raw URL (path is case-sensitive: IntPS not intps)
./tools/sd_to_pgbc_input.py \
  'https://raw.githubusercontent.com/symbolicdata/data/master/XMLResources/IntPS/Cyclic_4.xml' \
  -o input/from_sd_cyclic4.txt
```

Manifest: [`symbolicdata-manifest.yaml`](symbolicdata-manifest.yaml).

### Regression: converter vs committed inputs

For **cyclic4–7, cyclic9** and **katsura8**, the generator multiset in `input/*.txt` matches fresh conversion from SymbolicData `IntPS` (ordering of polynomials may differ; comparison sorts normalized strings). **Not** compared: `cyclic8` / `cyclic10` (no matching `Cyclic_3.xml` / `Cyclic_10.xml` in IntPS; in-tree `cyclic8` differs from `Cyclic_8.xml`), and **katsura7** (here: 7 polynomials; SD `Katsura_7.xml` is the 8-equation form in 8 variables).

```bash
make verify-symbolicdata   # needs network
```

### Adding the bundled “classic” benchmarks

Ten additional `input/*.txt` + `gb/*.txt` pairs (Caprasse, Butcher, …) are produced by:

```bash
make test
python3 tools/add_sd_benchmarks.py
```

Override binary or timeouts: `TEST_F4_BIN=...` `TIMEOUT=600`.

Run `test-f4` the same way as for other inputs:

```bash
./test/test-f4.bin input/from_sd_cyclic4.txt 1 0 1
```

## Caveats

- **`test-f4` uses GF(32003)**. Integer systems (`INTPS`) are reduced modulo 32003 when you run parallelGBC. For **`ModPS`**, if `<basedomain>` is not `GF(32003)`, the script still converts syntax but prints a **warning**: the numeric coefficients are meant for a different field.
- **Ordering**: SymbolicData does not always spell out the monomial order in these XML files. parallelGBC tests use **degrevlex** (`DegRevLexOrdering` in `test-f4`). Use the same order as in the literature for that benchmark when interpreting results.
- **Expected GB files** (`gb/*.txt`) are **not** produced by this script; generate them with your trusted CAS or `test-f4` once you trust the input.

## Supported XML shape

Root element **`INTPS`** or **`ModPS`** with:

- `<vars>a,b,c,...</vars>`
- `<basis><poly>...</poly>...</basis>`

Optional for `ModPS`: `<basedomain>GF(p)</basedomain>`.
