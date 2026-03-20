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
#pragma once
#include <vector>
#include "Polynomial.H"

/**
 * Print the given polynomials as coefficient matrix. The ordering defines
 * the column position of a given term.
 */
void printPolyMatrix(std::vector<Polynomial>& v, const TOrdering* O);
