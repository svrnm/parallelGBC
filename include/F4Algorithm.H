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
#include <unordered_map>
#include <algorithm>
#include <boost/bind.hpp>
#include "../include/F4Utils.H"

#ifdef F4_THREADS
	#define __F4_ALGORITHM_THREADS F4_THREADS
#else
	#define __F4_ALGORITHM_THREADS 1
#endif

class F4Pair {
	public:
		Term LCM;
		size_t i;
		size_t j;
		bool marked;
		degreeType sugar;

		F4Pair(const Term& LCM, size_t i, size_t j, bool marked, degreeType sugar) : LCM(LCM), i(i), j(j), marked(marked), sugar(sugar)	 {};
};

class F4PairComparator {
	public:
		const TOrdering* O;

		F4PairComparator(const TOrdering* O) : O(O) {}

		bool operator() (const F4Pair& lhs, const F4Pair& rhs) const 
		{
			return O->cmp(lhs.LCM, rhs.LCM) < 0;
		}
};

class F4Pivot {
	public:
		const Term column;
		size_t row;
		std::vector<std::pair<size_t, coeffType> > operations;

		F4Pivot(const Term& column, size_t row) : column(column), row(row) {}

		void addOperation(size_t row, coeffType factor) {
			operations.push_back(std::make_pair(row, factor));
		}
};

class F4Operation {
	public:
		size_t target;
		size_t oper;
		coeffType factor;

		F4Operation(size_t target, size_t oper, coeffType factor) : target(target), oper(oper), factor(factor) {}
};

typedef std::set<F4Pair,F4PairComparator> F4PairSet;

class F4 {
	public:
		CoeffField* field;
		std::vector<Polynomial> groebnerBasis;
		std::vector<bool> inGroebnerBasis;
		const TOrdering* O;
		degreeType currentDegree;
		int threads;

		// statistics
		double reductionTime;
		double prepareTime;
		double updateTime;

		// select() is part of reduce
		
		void gauss(std::vector<std::vector<coeffType> >& matrix, size_t upper, std::vector<bool>& empty);
		void pReduce(std::vector<std::vector<F4Operation> >& ops, std::vector<std::vector<coeffType> >& rs);
		size_t prepare(F4PairSet& pairs, std::vector<Polynomial>& polys, std::vector<std::vector<F4Operation> >& ops, std::set<Term, Term::comparator>& terms, std::vector<std::vector<coeffType> >& rs);
		void reduce(F4PairSet& pairs, std::vector<Polynomial>& polys);
		void postReduce(std::vector<Polynomial>& polys);
		void updatePairs(F4PairSet& pairs, std::vector<Polynomial>& polys, bool initial = false);
		std::vector<Polynomial> operator()(std::vector<Polynomial>& generators, const TOrdering* O, CoeffField* field, int threads = __F4_ALGORITHM_THREADS);
};

#endif
