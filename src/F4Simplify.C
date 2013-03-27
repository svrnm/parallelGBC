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

	

	std::pair<Term, Polynomial> F4Simplify::search(Term& t, Polynomial& f) {
		std::vector<degreeType> degree = t.getValues();
		std::vector<degreeType> helper(degree.size(), 0);
		bool check = true;
		while(check) {
			Term x(t.monoid(), helper);
			Term u = t.div ( x );
			for(size_t k = 0; k < F.size(); k++) {
				size_t j = F.size() - 1 - k;
				tbb::concurrent_unordered_map< Term, tbb::concurrent_vector<std::pair<Polynomial, Polynomial> > >::iterator it = F[j].find( u );
				if(it != F[j].end()) {
					for(size_t k = 0; k < it->second.size(); k++) {
						if(it->second[k].first == f) {
							Polynomial p = it->second[k].second;
							hits++;
							if( u.deg() == 0) {
								return std::make_pair( t, p );
							} else if( u.equal(t) ) {
								return std::make_pair( t.getOne(), p );
							} else {
								return search(x, p);
							}
						}
					}
				}
			}

			if(u.equal(t)) { 
				check = false; 
			} else {
				for(size_t i = 0; i < helper.size(); i++) {
					size_t j = helper.size() - 1 - i;
					helper[j]++;
					if(helper[j] > degree[j]) {
						helper[j] = 0;
					} else {
						break;
					}
				}
			}
		}
		misses++;
		return std::make_pair(t,f);
	}


	void F4Simplify::insert(size_t j, Term& t, Polynomial& f, Polynomial& p)  {
		if(f.size() >= p.size()) { return; }
		F[j][t].push_back(std::make_pair(f, p));
	}
}
