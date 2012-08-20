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
#include <set>
#include "../include/F4Utils.H"
#include "../include/Polynomial.H"

using namespace std;

void printPolyMatrix(vector<Polynomial>& v, const TOrdering* O)
{

	Term::comparator tog(O, true);
	set<Term, Term::comparator> terms(tog);


	for(size_t i = 0; i < v.size(); i++) {
		for(size_t j = 0; j < v[i].size(); j++) {
			terms.insert(v[i][j].second);
		}
	}


	for(size_t i = 0; i < v.size(); i++) {
		set<Term, Term::comparator>::iterator it = terms.begin();
		for(size_t j = 0; j < v[i].size(); j++) {
			for(;*it != v[i][j].second;it++) { std::cout << " 0 "; }
			std::cout << " " << v[i][j].first << " ";
			it++;
		}
		for(;it != terms.end(); it++) { std::cout << " 0 "; }
		std::cout << "\n";
	}

}

