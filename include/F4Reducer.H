 /**********************
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
			F4* f4;

			F4Reducer(F4* f4) : f4(f4) {}

			virtual ~F4Reducer() {}

			virtual void reduce(std::vector<Polynomial>& polys, degreeType currentDegree) = 0;
	
			virtual void addSPolynomial(size_t i, size_t j, Term& lcm) = 0;
	};
}
#endif
