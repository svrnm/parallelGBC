# parallelGBC

Parallel Gröbner basis computation

## About

This project implements a parallel algorithm for Gröbner basis computation. Background: [Groebner basis — Scholarpedia](https://www.scholarpedia.org/article/Groebner_basis).

The code comes from a master’s thesis and a paper in the *Proceedings of CASC 2012* (Maribor): [Springer chapter](http://link.springer.com/chapter/10.1007/978-3-642-32973-9_22).

## Requirements

- **CMake** 3.16+
- **C++20** compiler (e.g. GCC 10+, Clang 11+)
- **OpenMP** for C++ (e.g. GCC with `libgomp`). On **macOS + Apple Clang**, install **`libomp`** (`brew install libomp`); CMake uses `brew --prefix libomp` for flags. Without it, configure fails at `FindOpenMP`.
- [oneTBB](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onetbb.html) or compatible `libtbb`
- [Boost](https://www.boost.org/), especially **Boost.Regex** (for the `test/` driver)
- Multiple CPU cores if you want parallel speedups
- [SIMDe](https://github.com/simd-everywhere/simde) (git submodule) for portable SIMD (x86 / ARM)
- **Optional:** MPI + Boost.MPI for distributed runs — `-DENABLE_MPI=ON` (see [Build options](#build-options))

## Build and install

### Clone and submodules

```bash
git clone --recursive https://github.com/svrnm/parallelGBC
# already cloned:
git submodule update --init --recursive
```

### System packages (examples)

**Debian / Ubuntu** (similar to CI):

```bash
sudo apt-get install build-essential cmake g++ libboost-regex-dev libtbb-dev
```

**macOS** (Homebrew):

```bash
brew install cmake libomp boost tbb
```

### Configure, build, test

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

The example driver is **`build/test-f4`**.

### Install (optional)

```bash
cmake --install build --prefix /usr/local
```

### Build options

Defaults match the historical in-tree defaults. Set at configure time, e.g. `cmake -B build -DPGBC_USE_SSE=0`, or use `ccmake` / `cmake-gui`:

| Variable | Meaning |
| --- | --- |
| `PGBC_COEFF_BITS` | Coefficient storage bits (8, 16, 32) |
| `PGBC_USE_SSE` | SIMDe SIMD field ops (`0` or `1`) |
| `PGBC_SORTING` | Reduction sorting mode |
| `PGBC_POST_REDUCE` | Post-reduce with simplify |
| `PGBC_PARALLEL_SETUP` | Parallel matrix setup |
| `ENABLE_MPI` | MPI + Boost.MPI (`ON` / `OFF`, default `OFF`; requires MPI C++ and Boost.MPI/serialization) |

### Using the installed CMake package

After install: `find_package(parallelGBC CONFIG)` and `target_link_libraries(... parallelGBC::f4)`. Headers are under **`include/parallelGBC/`** (e.g. `#include <parallelGBC/F4.H>`).

### Continuous integration

[`.github/workflows/ci.yml`](.github/workflows/ci.yml) configures with CMake, runs **`ctest`**, and smoke-tests install (**`<prefix>/lib/libf4.a`** after `cmake --install`).

## Testing

Example: degree reverse lexicographic Gröbner basis of **cyclic-8**, 4 threads, high verbosity, no GB printout, block size 1024, no simplify, with sugar — from the repository root:

```bash
./build/test-f4 input/cyclic8.txt 4 127 0 1024 0 1
```

General invocation:

```text
./build/test-f4 <input-file> <processors> <verbosity> <printGB> <blocksize> <doSimplify> <withSugar>
```

With MPI:

```bash
mpirun -np <slots> --host <hosts> ./build/test-f4 <...>
```

## Regression / reference outputs

The **`gb/`** directory holds reference Gröbner bases over **F₃₂₀₀₃** in degree reverse lex order (ApCoCoA-style).

Run the full regression suite:

```bash
ctest --test-dir build --output-on-failure
```

Or the shell driver (from repo root; override binary if needed):

```bash
TEST_F4_BIN="$PWD/build/test-f4" ./test/RunTests.sh
```

Behavior:

- Compares against **`gb/*.txt`** on the **smallest** value in `CORE_LIST`.
- On the **largest** core count: parallel F4 and, by default, Buchberger verification (OpenMP S-pairs), unless expected |G| exceeds `VERIFY_MAX_GB`.

Useful environment overrides:

| Scenario | Example |
| --- | --- |
| Skip Buchberger verification | `VERIFY_GB=0 ./test/RunTests.sh` |
| Verify all sizes (can be very slow) | `VERIFY_MAX_GB=0` |
| Cap verification time (seconds) | `VERIFY_TIMEOUT=900` |

Defaults in **`test/RunTests.sh`**: `CORE_LIST` is `1 8`; `VERIFY_GB` is `0` if unset (use `VERIFY_GB=1` for stricter runs). **`ctest`** / CI set `VERIFY_GB`, `VERIFY_MAX_GB`, `VERIFY_TIMEOUT`, and `CORE_LIST` for you.

## Developer tooling

- **`.clang-format`** — optional for `.C` / `.H`; run from your editor or CLI.
- **`compile_commands.json`** — for clangd / IDEs: configure with `cmake -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON`, build once, then point the IDE at **`build/compile_commands.json`**, or run **`./scripts/gen-compile_commands.sh`** to symlink it to the repo root (gitignored).

## Verbosity flags

Runtime verbosity does not materially affect performance. It is an extra argument to the F4 operator; you can also pass an output stream. Default: no verbosity, `std::cout`.

| Bit | Meaning |
| --- | --- |
| 1 | Runtime |
| 2 | Reduction time |
| 4 | Prepare time |
| 8 | Update time |
| 16 | Print sugar degree during reduction |
| 32 | Print time of reduction step |
| 64 | Print matrix size during reduction |
| 128 | Print all computed polynomials |

## Usage in your own project

See [`test/test-f4.C`](test/test-f4.C) for a minimal example of linking against the library.

## Contact

Severin Neumann — [severin.neumann@altmuehlnet.de](mailto:severin.neumann@altmuehlnet.de)

## License

This program is free software; see [LICENSE.txt](LICENSE.txt) for details.
