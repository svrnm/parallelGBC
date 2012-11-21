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

	// Setup a term comparator which return true if a term a is greater than a term b
	Term::comparator tog(O, true);

	// Set which collects all terms used in v.
	set<Term, Term::comparator> terms(tog);


	// Read in all terms.
	for(size_t i = 0; i < v.size(); i++) {
		for(size_t j = 0; j < v[i].size(); j++) {
			terms.insert(v[i][j].second);
		}
	}

	// Print the matrix.
	for(size_t i = 0; i < v.size(); i++) {
		set<Term, Term::comparator>::iterator it = terms.begin();
		// Iterate over all values in the current row
		for(size_t j = 0; j < v[i].size(); j++) {
			// Each term which does not occure in the current polynomial
			// is zero.
			for(;*it != v[i][j].second;it++) { std::cout << " 0 "; }
			// Output the current value
			std::cout << " " << v[i][j].first << " ";
			it++;
		}
		// If there are more terms which are less than the smallest term
		// within the current polynomial, print zeros for each of these terms.
		for(;it != terms.end(); it++) { std::cout << " 0 "; }

		// Break line after one polynomial has been processed
		std::cout << "\n";
	}
}

