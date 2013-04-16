 /**
 *  An abstraction of the 'reduction' part of the F4 algorithm. This allows to plug-in further
 *  reducers.
 *
 ****
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
#ifndef F4_REDUCER_H
#define F4_REDUCER_H
#include <set>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/bind.hpp>
#include <tbb/concurrent_unordered_map.h>
#include <unordered_map>
#include <tbb/concurrent_vector.h>
#include <tbb/blocked_range.h>
#include <tbb/atomic.h>
#include "../include/Term.H"
#include "../include/Polynomial.H"
#include "../include/F4Logger.H"

namespace parallelGBC {

	class F4;

	class F4Reducer {
		public:
			// The reducer is directly linked to one instance of 'F4'. So you have access to all
			// the other properties (groebnerBasis, inGroebnerBasis, logger, ...)
			F4* f4;

			// Construct a new instance of F4Reducer
			F4Reducer(F4* f4) : f4(f4) {}

			virtual ~F4Reducer() {}

			// Overwrite if required. This method is called before the first reduction
			virtual void init() { }

			// MPI specific (?) for non-master reducers, implement if required
			virtual void client() { }

			// Overwirte if required. This method is called after the last reduction, e.g. for clean up
			virtual void finish() { }

			// Implement to realize your specific reduction
			virtual void reduce(std::vector<Polynomial>& polys, degreeType currentDegree) = 0;
	
			// This function is called to provide the list of all S-polynomials. In one reduction there may
			// be several s-polynomials to be reduced
			virtual void addSPolynomial(size_t i, size_t j, Term& lcm) = 0;
	};
}
#endif
