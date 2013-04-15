/**
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
#ifndef F4_SimplifyDB_H
#define F4_SimplifyDB_H
#include <vector>
#include <unordered_map>
#include <tbb/concurrent_unordered_map.h>
#include "../include/Term.H"
#include "../include/Polynomial.H"
#include "../include/TOrdering.H"

namespace parallelGBC {

	class F4SimplifyDB {
		protected:
			const TOrdering* O;
			tbb::concurrent_unordered_map<size_t, tbb::concurrent_unordered_map<Term, Polynomial, std::hash<Term> > > database;

		public:
			F4SimplifyDB(const TOrdering* O) : O(O) { }

			std::pair<Term, Polynomial> search(size_t i, Term& t, bool full = false);

			size_t check(size_t i, Term& t);
			
			void insert(size_t i, Term& t, Polynomial& p);
			
			~F4SimplifyDB() {
				database.clear();
			}

	};

}
#endif
