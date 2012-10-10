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
#include <algorithm>
#include <boost/bind.hpp>
#include <iostream>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_vector.h>
#include <tbb/blocked_range.h>
#include "../include/Term.H"
#include "../include/Polynomial.H"
#include "../include/F4Utils.H"

namespace parallelGBC {

	struct F4Pair {
		Term LCM;
		size_t i;
		size_t j;
		bool marked;
		degreeType sugar;

		F4Pair(const Term& LCM, size_t i, size_t j, bool marked, degreeType sugar) : LCM(LCM), i(i), j(j), marked(marked), sugar(sugar)	 {};

		struct comparator {
			const TOrdering* O;

			comparator() : O(NULL) {}
			comparator(const TOrdering* O) : O(O) {}

			bool operator() (const F4Pair& lhs, const F4Pair& rhs) const 
			{
				return O->cmp(lhs.LCM, rhs.LCM) < 0;
			}
		};

		struct sugarComparator
		{
			bool operator() (const F4Pair& a, const F4Pair& b)
			{
				return a.sugar < b.sugar;
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
		std::vector<int32_t> targets;
		std::vector<int32_t> opers;
		coeffRow factors;

		void push_back(int32_t t, int32_t o, coeffType f) {
			targets.push_back( t );
			opers.push_back( o );
			factors.push_back( f );
		}

		size_t size() {
			return targets.size();
		}
	};


	typedef std::multiset<F4Pair,F4Pair::comparator> F4PairSet;


	class F4 {
		public:


			std::vector<Polynomial> groebnerBasis;
			std::vector<bool> inGroebnerBasis;
			const TOrdering* O;
			CoeffField* field;
			int threads;
			degreeType currentDegree;

			// Datastructures for computation:
			F4PairSet pairs;
			std::set<Term, Term::comparator> terms;

			tbb::concurrent_vector<std::pair<size_t, Term> > rows;
			tbb::concurrent_unordered_map<Term, uint32_t, std::hash<Term> > pivots;
			std::map<Term, uint32_t, Term::comparator> pivotsOrdered;
			tbb::concurrent_unordered_map<Term, std::vector<std::pair<uint32_t, coeffType> >, std::hash<Term> > pivotOps;
			tbb::concurrent_unordered_set<Term, std::hash<Term> > termsUnordered;
			tbb::concurrent_vector<tbb::concurrent_vector<Monomial> > rightSide;
			std::vector<F4Operations > ops;
			coeffMatrix rs;
			std::vector<size_t> deps;

			std::vector<size_t> prefixes;
			std::vector<size_t> suffixes;

			size_t upper;

			std::vector<bool> empty;

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

			// Helper functions for groebner bases computation
			void gauss();
			void pReduce();
			void prepare(std::vector<Polynomial>& polys);
			void reduce(std::vector<Polynomial>& polys);
			void updatePairs(std::vector<Polynomial>& polys);

			void setupRow(Polynomial& current, Term& ir, size_t i, tbb::blocked_range<size_t>& range);
			void setupDenseRow(tbb::blocked_range<size_t>& range);
			void pReduceRange(size_t i, tbb::blocked_range<size_t>& range);


			// Call this operator to compute the groebner basis of "generators"
			std::vector<Polynomial> operator()(std::vector<Polynomial>& generators, const TOrdering* O, CoeffField* field, int threads = 1, int verbosity = 0, std::ostream& output = std::cout);
	};

	struct F4SetupRow
	{
		F4& f4;
		Polynomial& current;
		Term& ir;
		size_t i;

		F4SetupRow(F4& f4, Polynomial& current, Term& ir, size_t i) : f4(f4), current(current), ir(ir), i(i) {}

		void operator() (tbb::blocked_range<size_t>& range) const
		{
			f4.setupRow(current, ir, i, range);
		}
	};

	struct F4SetupDenseRow {
		F4& f4;

		F4SetupDenseRow(F4& f4) : f4(f4) {}

		void operator() (tbb::blocked_range<size_t>& range) const
		{
			f4.setupDenseRow(range);
		}
	};

	struct F4PReduceRange
	{
		F4& f4;
		size_t i;

		F4PReduceRange(F4& f4, size_t i) : f4(f4), i(i) {}

		void operator() (tbb::blocked_range<size_t>& range) const
		{
			f4.pReduceRange(i, range);
		}
	};



}
#endif
