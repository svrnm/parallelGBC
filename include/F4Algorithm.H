#ifndef F4_ALGORITHM_H
#define F4_ALGORITHM_H
#include <set>
#include <vector>
#include <map>
#include "../include/Term.H"
#include "../include/Polynomial.H"
#include <boost/unordered_map.hpp>
#include <boost/bind.hpp>
#include "../include/F4Utils.H"

#ifdef F4_THREADS
	#define __F4_ALGORITHM_THREADS F4_THREADS
#else
	#define __F4_ALGORITHM_THREADS 1
#endif



using namespace boost;

class F4Pair {
	public:
		const Term* LCM;
		size_t i;
		size_t j;
		bool marked;
		degreeType sugar;

		F4Pair(const Term* LCM, size_t i, size_t j, bool marked, degreeType sugar) : LCM(LCM), i(i), j(j), marked(marked), sugar(sugar)	 {};
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
		const Term* column;
		size_t row;
		vector<pair<size_t, coeffType> > operations;

		F4Pivot(const Term* column, size_t row) : column(column), row(row) {}

		void addOperation(size_t row, coeffType factor) {
			operations.push_back(make_pair(row, factor));
		}
};

class F4Operation {
	public:
		size_t target;
		size_t oper;
		coeffType factor;

		F4Operation(size_t target, size_t oper, coeffType factor) : target(target), oper(oper), factor(factor) {}
};

typedef set<F4Pair,F4PairComparator> F4PairSet;

class F4 {
	public:
		CoeffField* field;
		vector<Polynomial> groebnerBasis;
		vector<bool> inGroebnerBasis;
		const TOrdering* O;
		degreeType currentDegree;
		int threads;

		// statistics
		double reductionTime;
		double prepareTime;
		double updateTime;

		// select() is part of reduce
		
		void gauss(vector<vector<coeffType> >& matrix, size_t upper, vector<bool>& empty);
		void pReduce(vector<vector<F4Operation> >& ops, vector<vector<coeffType> >& rs);
		size_t prepare(F4PairSet& pairs, vector<Polynomial>& polys, vector<vector<F4Operation> >& ops, set<const Term*, TermComparator>& terms, vector<vector<coeffType> >& rs);
		void reduce(F4PairSet& pairs, vector<Polynomial>& polys);
		void postReduce(vector<Polynomial>& polys);
		void updatePairs(F4PairSet& pairs, vector<Polynomial>& polys, bool initial = false);
		vector<Polynomial> operator()(vector<Polynomial>& generators, const TOrdering* O, CoeffField* field, int threads = __F4_ALGORITHM_THREADS);
};

#endif
