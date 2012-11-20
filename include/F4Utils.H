/**
 *  Provides helper functions which are used for timings and debugging:
 *
 *  - seconds() can be used to measure the time of a function:
 *		double s = seconds();
 *		doSomething();
 *		double time = seconds()-s;
 *
 *	- printPolyMatrix(...) writes a vector of polynomials as coefficient matrix to
 *	the standard output.
 *
 ***************
 *
 *  This file is part of parallelGBC, a parallel groebner basis computation tool.
 *
 *  parallelGBC is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  parallelGBC is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with parallelGBC.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef F4_UTILS_H
#define F4_UTILS_H
#include <sys/time.h>
#include <iostream>
#include <vector>
#include <sys/types.h>
#include <sys/sysinfo.h>
#include "Polynomial.H"

/**
 * Return the current time stamp in seconds.
 */
inline double seconds()
{    
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + ((double)tv.tv_usec)/1000000.0;
}

/**
 * Print the given polynomials as coefficient matrix. The ordering defines
 * the column position of a given term.
 */
void printPolyMatrix(std::vector<Polynomial>& v, const TOrdering* O);
#endif
