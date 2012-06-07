Parallel Groebner Basis Computation (alpha 0.1)

LICENSE
This program is free software; see LICENSE.TXT for more details

REASON
This program provides an algorithm for parallel groebner basis computation.
If you do not know, what a groebner basis is and what they are for,
you may read on here:

http://www.scholarpedia.org/article/Groebner_basis

The code of the project is the result of my master thesis. Documents for
reading will follow in later releases (some translations have to be done...)


REQUIREMENTS
* A compiler which supports OpenMP and C++0x / C++11 (GCC is fine)
* A processor which has SSE2, if not disable it in the Makefile.rules
* Several processors if you want to use the parallelization (dual or quadcores, etc.)

INSTALLATION
Just execute

# make

and if you have several CPUs, you can use

# make -j<NUM_OF_PROCS>

If you experience any problems, then have a look to the Makefile or 
contact the author

TESTING
Compute the degree reverse lexicogrpahic gröbner basis of cyclic-8 with 4 threads
and a lot of verbosity

cd test/;
./test/test-f4 ../input/cyclic8.txt 4 127

VERBOSITY
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

USAGE AND EXAMPLE
If you want to use the code for your own project see test/test-f4.C as example

CONTACT
For any questions you can contact Severin Neumann <severin.neumann@computer.org>
