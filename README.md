Parallel Groebner Basis Computation (beta 0.9)
=======

License
-------
This program is free software; see LICENSE.txt for more details

Reason
------
This program provides an algorithm for parallel groebner basis computation.
If you do not know, what a groebner basis is and what they are for,
you may read on here [Scholarpedia](https://www.scholarpedia.org/article/Groebner_basis).

The code of the project is the result of my master thesis and a paper published to the
Proceedings of CASC 2012 in Maribor. You can read the paper at [Springer Link](http://link.springer.com/chapter/10.1007/978-3-642-32973-9_22).


Requirements
------------
* A compiler with **C++20** support (e.g. GCC 10+, Clang 11+)
* [oneAPI Threading Building Blocks (oneTBB)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onetbb.html) or compatible `libtbb`
* [Boost](http://www.boost.org/), especially Boost.Regex (if you want to use the example binaries in test/)
* OpenMP is optional but can speed up some more computations by parallelization
* Several processors if you want to use the parallelization (dual or quadcores, etc.).
* [SIMDe](https://github.com/simd-everywhere/simde) (included as git submodule) for portable SIMD on x86 and ARM.
* openmpi and Boost.MPI if you want to do distributed parallelization, if not disable the MPI option in Makefile.rules.

Installation
------------
Clone the repository and initialize the SIMDe submodule:

    git clone --recursive <repository-url>
    # or, if already cloned:
    git submodule update --init

If you need to configure some settings (SSE, MPI) just have a look into Makefile.rules

		vim Makefile.rules

Afterwards or if you'd like to use the default settings just execute

    make

and if you have several CPUs, you can use

    make -j<NUM_OF_PROCS>

If you experience any problems, then have a look to the Makefile.rules or contact the author

Testing
-------
Compute the degree reverse lexicographic gröbner basis of cyclic-8 with 4 threads
and a lot of verbosity and without printing the groebner basis. The block size of
the matrix is 1024. For computation the simplify algorithm is not used, the sugar
cube selection strategy is.

    ./test/test-f4.bin ../input/cyclic8.txt 4 127 0 1024 0 1

In general you can compute with this binary using the following parameters:

    ./test/test-f4.bin <input-file> <processors> <verbosity> <printGB> <blocksize> <doSimplify> <withSugar>

If you have compiled the binary using MPI you can compute distributed:

		mpirun -np <slots> --host <hosts> ./test/test-f4.bin <...>

Checking functionality
----------------------
The folder gb/ contains precomputed groebner bases over F_{32003} using degree reverse
lexicographic term ordering (computed using ApCoCoA). Use 

    'make check'
        
to validate the functionality of parallelGBC.

`make check` runs Buchberger verification by default on the **largest** thread count in `CORE_LIST` (`1 8` unless you override), only when the expected basis size |G| is at most `VERIFY_MAX_GB` (default **100**). S-pair reduction uses **OpenMP** inside the verifier. Use `VERIFY_GB=0 make check` to skip verification. Set `VERIFY_MAX_GB=0` to verify all sizes (can be very slow). Each verify run is capped by `VERIFY_TIMEOUT` seconds (default **600** for `make check`).

Verbosity
---------
Verbosity, which can be changed during runtime, nothing which should
influence performance. It' is an additional parameter for the F4 operator().
Additionally you can give an output stream, which should be used for output.
Default ist no verbosity and std::cout as output stream.

1 - Runtime

2 - Reduction time

4 - Prepare time

8 - Update time

16 - Print sugar degree during reduction step

32 - Print time of reduction step

64 - Print matrix size during reduction step

128 - Print all computed polynomials

Usage and example
-----------------
If you want to use the code for your own project see test/test-f4.C as example.

Contact
-------
For any questions you can contact Severin Neumann <severin.neumann@computer.org>
