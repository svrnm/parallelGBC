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
#include "../include/F4Simplify.H"
namespace parallelGBC {



	//	std::pair<Term, Polynomial> 
	void F4Simplify::search(Term& t, Polynomial& f) {
		tbb::concurrent_unordered_map<Polynomial, tbb::concurrent_unordered_map<Term, Polynomial, std::hash<Term> > >::iterator it = F.find( f );
		if(it != F.end()) {
			std::vector<degreeType> degree = t.getValues();
			std::vector<degreeType> helper(degree.size(), 0);
			bool check = true;
			while(check) {
				Term x(t.monoid(), helper);
				Term u = t.div ( x );
				tbb::concurrent_unordered_map<Term, Polynomial>::iterator it2 = it->second.find( u );
				if(it2 != it->second.end()) {
					f = it2->second;
					if( u.deg() == 0) {
						return;
					} else if( u.equal(t) ) {
						t = t.getOne();
						return;
					} else {
						t = x;
						return search(t, f);
					}
				}

				if(u.equal(t)) { 
					check = false; 
				} else {
					for(size_t j = 0; j < helper.size(); j++) {
						helper[j]++;
						if(helper[j] > degree[j]) {
							helper[j] = 0;
						} else {
							break;
						}
					}
				}
			}
		}
		//return std::make_pair(t,f);
	}


	void F4Simplify::insert(Term& t, Polynomial& f, Polynomial& p)  {
		if(f == p) { return; }
		F[f][t] = p;
	}
}
