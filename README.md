Parallel Groebner Basis Computation (beta 0.8)
=======

License
-------
This program is free software; see LICENSE.TXT for more details

REASON
------
This program provides an algorithm for parallel groebner basis computation.
If you do not know, what a groebner basis is and what they are for,
you may read on here [Scholarpedia](http://www.scholarpedia.org/article/Groebner_basis).

The code of the project is the result of my master thesis and a paper published to the
Proceedings of CASC 2012 in Maribor. You can read the paper at [Springer Link](http://link.springer.com/chapter/10.1007/978-3-642-32973-9_22).


Requirements
------------
* A compiler which supports C++0x / C++11 (GCC is fine)
* [Intel TBB](http://threadingbuildingblocks.org/)
* A processor which has SSE2, if not disable SSE in the Makefile.rules
* Several processors if you want to use the parallelization (dual or quadcores, etc.)

Installation
------------
Just execute

    make

and if you have several CPUs, you can use

    make -j<NUM_OF_PROCS>

If you experience any problems, then have a look to the Makefile or 
contact the author

Testing
-------
Compute the degree reverse lexicogrpahic gröbner basis of cyclic-8 with 4 threads
and a lot of verbosity and without printing the groebner basis. The block size of
the matrix is 1024. For computation the simplify algorithm is not used, the sugar
cube selection strategy is.

    ./test/test-f4 ../input/cyclic8.txt 4 127 0 1024 0 1

In general you can compute with this binary using the following parameters:

    ./test/test-f4 <input-file> <processors> <verbosity> <printGB> <blocksize> <doSimplify> <withSugar>

Checking functionality
----------------------
The folder gb/ contains precomputed groebner bases over F_{32003} using degree reverse
lexicographic term ordering (computed using ApCoCoA). Use 

    'make check'
        
to validate the functionality of parallelGBC.

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
