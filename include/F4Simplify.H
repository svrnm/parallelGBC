/**
 **********************
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
#ifndef F4_SIMPLIFY_H
#define F4_SIMPLIFY_H
#include <sys/time.h>
#include <iostream>
#include <sys/types.h>
#include <sys/sysinfo.h>
#include <vector>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_vector.h>
#include "../include/Polynomial.H"

namespace parallelGBC {

	class F4Simplify {


		public:
			tbb::concurrent_unordered_map<Polynomial, tbb::concurrent_unordered_map<Term, Polynomial, std::hash<Term> >, std::hash<Polynomial> > F;

			F4Simplify() { }

			void search(Term& t, Polynomial& f);

			void insert(Term& t, Polynomial& f, Polynomial& p);
	};
}
#endif
