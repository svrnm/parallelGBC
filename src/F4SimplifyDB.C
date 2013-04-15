/*
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
#include "../include/F4SimplifyDB.H"

using namespace std;

namespace parallelGBC {

	std::pair<Term, Polynomial> F4SimplifyDB::search(size_t i, Term& t, bool full) {
		Term u = t;
		// Best case is if there is alreade a t*f stored.
		tbb::concurrent_unordered_map<Term, Polynomial>::iterator p = database[i].find(t);

		if(p == database[i].end()) {
			if(full) {
				for(tbb::concurrent_unordered_map<Term, Polynomial>::iterator it = database[i].begin(); it != database[i].end(); it++) {
					if(t.isDivisibleBy(it->first)) {
						Term x = t.div( it->first ); // x = t/u;
						if( O->cmp(u, x) >= 0 ) { // >= takes also care of p: p should not be empty!
							u = x;
							p = it;
						}
					}
				}
			} else {
				vector<Term> divisors = t.divAllX();
				for(size_t j = 0; j < divisors.size() && p == database[i].end(); j++) {
					p = database[i].find(divisors[j]);
				}
				if(p == database[i].end()) {
					p = database[i].find(t.getOne());
					u = t;
				} else {
					u = t.div( p->first );
				}
			}
		} else {
			u = t.getOne();
		}
		return std::make_pair(u, p->second);
	}

	size_t F4SimplifyDB::check(size_t i, Term& t) {
		tbb::concurrent_unordered_map<Term, Polynomial>::iterator it = database[i].find( t );
		if( it == database[i].end() ) {
			return (size_t)-1;
		} else {
			return database[i][t].size();
		}
	}

	void F4SimplifyDB::insert(size_t i, Term& t, Polynomial& p) {
		tbb::concurrent_unordered_map<Term, Polynomial>::iterator it = database[i].find( t );
		if( it == database[i].end() ) {
			database[i].insert( std::make_pair(t,p) );
		} else  {
			database[i][t] = p;
		}
	}
}
