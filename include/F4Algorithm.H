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
#ifndef F4_ALGORITHM_H
#define F4_ALGORITHM_H
#include <set>
#include <vector>
#include <map>
#include "../include/Term.H"
#include "../include/Polynomial.H"
#include <algorithm>
#include <boost/bind.hpp>
#include <iostream>
#include "../include/F4Utils.H"

struct F4Pair {
		Term LCM;
		size_t i;
		size_t j;
		bool marked;
		degreeType sugar;

		F4Pair(const Term& LCM, size_t i, size_t j, bool marked, degreeType sugar) : LCM(LCM), i(i), j(j), marked(marked), sugar(sugar)	 {};

	struct comparator {
		const TOrdering* O;

		comparator(const TOrdering* O) : O(O) {}

		bool operator() (const F4Pair& lhs, const F4Pair& rhs) const 
		{
			return O->cmp(lhs.LCM, rhs.LCM) < 0;
		}
	};
};

struct F4Pivot {
		const Term column;
		size_t row;
		std::vector<std::pair<size_t, coeffType> > operations;

		F4Pivot(const Term& column, size_t row) : column(column), row(row) {}

		void addOperation(size_t row, coeffType factor) {
			operations.push_back(std::make_pair(row, factor));
		}
};

struct F4Operations {
		std::vector<size_t> targets;
		std::vector<size_t> opers;
		std::vector<coeffType> factors;

		void push_back(size_t t, size_t o, coeffType f) {
			targets.push_back( t );
			opers.push_back( o );
			factors.push_back( f );
		}

		size_t size() {
			return targets.size();
		}
};

typedef std::set<F4Pair,F4Pair::comparator> F4PairSet;

class F4 {
	public:
		std::vector<Polynomial> groebnerBasis;
		std::vector<bool> inGroebnerBasis;
		const TOrdering* O;
		CoeffField* field;
		int threads;
		degreeType currentDegree;

		/* Verbosity, which can be changed during runtime, nothing which should
		 * influence performance.
		 *
		 * 1 - Runtime
		 * 2 - Reduction time
		 * 4 - Prepare time
		 * 8 - Update time
		 * 16 - Print sugar degree during reduction step
		 * 32 - Print time of reduction step
		 * 64 - Print matrix size during reduction step 
		 */
		int verbosity;
		// There we write verbose messages
		std::ostream* out;

		// statistics
		double reductionTime;
		double prepareTime;
		double updateTime;

		// select() is part of reduce
		
		void gauss(std::vector<std::vector<coeffType> >& matrix, size_t upper, std::vector<bool>& empty);
		void pReduce(std::vector<F4Operations>& ops, std::vector<size_t>& deps, std::vector<std::vector<coeffType> >& rs);
		size_t prepare(F4PairSet& pairs, std::vector<Polynomial>& polys, std::vector<F4Operations>& ops, std::vector<size_t>& deps, std::set<Term, Term::comparator>& terms, std::vector<std::vector<coeffType> >& rs);
		void reduce(F4PairSet& pairs, std::vector<Polynomial>& polys);
		void updatePairs(F4PairSet& pairs, std::vector<Polynomial>& polys, bool initial = false);
		std::vector<Polynomial> operator()(std::vector<Polynomial>& generators, const TOrdering* O, CoeffField* field, int threads = 1, int verbosity = 0, std::ostream& output = std::cout);
};

#endif
