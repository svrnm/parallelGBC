/**
 * This file includes the headers for the parallel implementation of the f4 algorithm which
 * computes groebner basis for a given term ordering. The class "F4" provides a operator to
 * execute the algorithm for a given set of polynomials, a term ordering and a coefficient field 
 * using a given number of parallel threads. A usage example is given in test/test-f4.C and in
 * test/test-lex.C
 *
 **********************
 *
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
#include <tbb/atomic.h>
#include "../include/Term.H"
#include "../include/Polynomial.H"
#include "../include/F4Utils.H"

namespace parallelGBC {

	/**
	 * This class represents a critical pair during the groebner basis computation.
	 * A critical pair is identified by the indices of the elements of the already
	 * found groebner basis elements. Further more the least common multiple of the
	 * leading terms of these two polynomials and the sugar degree of the resulting
	 * S-polynomial is cached. The flag 'marked' is used during the UPDATE(...) function
	 * as described by Gebauer and Möller (http://dl.acm.org/citation.cfm?id=65042)
	 */
	struct F4Pair {
		/**
		 * The least common multiple of LT(f_i) and LT(f_j)
		 */
		Term LCM;
		/**
		 * The index of the first groebner basis element
		 */
		size_t i;
		/**
		 * The index of the second groebner basis element
		 */
		size_t j;
		/**
		 * This flag is set to true if LCM == groebnerBasis[i].LT().mul(h.LT()) which satisfies the third criterium of the
		 * Gebauer Möller installation ("...Finally delete all (i,t) with T(i)T(j) = T(i,t)...")
		 */
		bool marked;
		/**
		 * The sugar of the S-polynomial 
		 */
		degreeType sugar;

		/**
		 * Construct an instance of F4Pair. All fields are set at setup time
		 */
		F4Pair(const Term& LCM, size_t i, size_t j, bool marked, degreeType sugar) : LCM(LCM), i(i), j(j), marked(marked), sugar(sugar)	 {};

		/**
		 * This class represents a comparator for F4Pairs. An instance of F4Pair is less than another F4Pair if the LCM of the
		 * first pair is less than the second. To check if 'less than' a term ordering O is given.
		 */
		struct comparator {
			/**
			 * The term ordering which is used to compare the pairs
			 */
			const TOrdering* O;

			/**
			 * Copy constructor. Attention, this may lead to misbehaviour if
			 * O is not set to a real term ordering.
			 */
			comparator() : O(NULL) {}

			/**
			 * Construct an instance of F4Pair::comparator using the given term ordering O
			 */
			comparator(const TOrdering* O) : O(O) {}

			/**
			 * Returns true if the LCM of lhs is less than the LCM of rhs
			 */
			bool operator() (const F4Pair& lhs, const F4Pair& rhs) const 
			{
				return O->cmp(lhs.LCM, rhs.LCM) < 0;
			}
		};

		/**
		 * This class represents a comparator for F4Pairs which depends on the sugar degree of two tiven F4Pairs
		 * A F4Pair a is less than a F4Pair b if the sugar degree of a is less than the sugar degree of b.
		 */
		struct sugarComparator
		{
			/**
			 * Return true if the sugar degree of a is less than the sugar degree of b
			 */
			bool operator() (const F4Pair& a, const F4Pair& b)
			{
				return a.sugar < b.sugar;
			}
		};


	};

	/**
	 * This class represents a container for matrix operations which can be computed in parallel. An Operation
	 * is a triple of a target row 't', an operator row 'o' and a coefficient 'f', which means that 'o' times 'f'
	 * should be subtracted from 't'. 
	 */
	struct F4Operations {
		/**
		 * The list of targets. A target at position i belongs to a operator and a factor at position i
		 */
		std::vector<int32_t> targets;
		/**
		 * The list of operators. A operator at position i belongs to a target and a factor at position i.
		 */
		std::vector<int32_t> opers;
		/**
		 * The list of factors. A factor at position i belongs to a operator and a coefficient at position i.
		 */
		coeffRow factors;

		/**
		 * Add a new operation (target,operator,factor) to the given container of operations. 
		 */
		void push_back(int32_t t, int32_t o, coeffType f) {
			targets.push_back( t );
			opers.push_back( o );
			factors.push_back( f );
		}

		int32_t target(size_t i) {
			return targets[i];
		}

		int32_t oper(size_t i) {
			return opers[i];
		}

		coeffType factor(size_t i) {
			return factors[i];
		}

		/**
		 * Returns the size of the container which is equivalent to the number of stored operations
		 */
		size_t size() {
			return targets.size();
		}
	};

	/**
	 * Shortname for a multiset containing F4Pairs. A multiset is used, since it possible that
	 * two pairs have the same least common multiple.
	 */
	typedef std::multiset<F4Pair,F4Pair::comparator> F4PairSet;


	/**
	 * This class represents the F4 algorithm. It is implemented as function object, therefor you
	 * can create a instance of F4 'f4' and call the algorithm by 'result = f4(generators,O,field,threads,verbosity,output);'.
	 * The parameters are:
	 *
	 * - generators: The input system of polynomials for which a groebner basis should be computed
	 * - O: The term ordering which should be used to compute the groebner basis
	 * - field: The coefficient field for the computations during reduction.
	 * - threads: The number of parallel threads used to compute the solution
	 * - verbosity: The level of debugging informations which are printed to 'output'
	 * - output: The stream which writes the debug informations out.
	 */
	class F4 {
		public:
			/**
			 * The (intermediate) result
			 */
			std::vector<Polynomial> groebnerBasis;
			/**
			 * For each element in 'groebnerBasis' the information is stored, if the
			 * element i will be part of the final result
			 */
			std::vector<bool> inGroebnerBasis;
			/**
			 * The used term ordering
			 */
			const TOrdering* O;
			/**
			 * The used coefficient fiel
			 */
			CoeffField* field;
			/**
			 * The number of computation thread
			 */
			int threads;

			/**
			 * The current sugar degree, which will (or has been) computed
			 */
			degreeType currentDegree;

			/**
			 * The set of criticial pairs which will be computed
			 */
			F4PairSet pairs;
			/**
			 * The set of terms which occure during reduction and are not leading terms
			 */
			std::map<Term, uint32_t, Term::comparator> terms;

			/**
			 * Intermediate storage for critical pairs which will be used for the next
			 * reduction. For each critical pair there will be two elements in the list
			 * rows identified by their index of the existing groebner basis and the lcm
			 * of the two pairs. 
			 */
			tbb::concurrent_vector<std::pair<size_t, Term> > rows;
			/**
			 * The set of terms which occure during reduction and are leading terms (=pivots).
			 * The value of the map is the row which is the pivot row.
			 */
			tbb::concurrent_unordered_map<Term, uint32_t, std::hash<Term> > pivots;
			/**
			 * A ordered copy of 'pivots'
			 */
			std::map<Term, uint32_t, Term::comparator> pivotsOrdered;
			/**
			 * Intermediate storage of the operations which have to be executed to reduce the matrix
			 */
			tbb::concurrent_unordered_map<Term, std::vector<std::pair<uint32_t, coeffType> >, std::hash<Term> > pivotOps;
			/**
			 * A unordered copy of 'terms'
			 */
			tbb::atomic<uint32_t> termCounter;
			tbb::concurrent_unordered_map<Term, uint32_t, std::hash<Term> > termsUnordered;
			
			/**
			 * A sparse representation of the matrix entries in the non pivot part
			 */
			//tbb::concurrent_unordered_map<uint32_t, tbb::concurrent_vector<std::pair<coeffType, uint32_t> > > rightSide;
			tbb::concurrent_vector<tbb::concurrent_vector<std::pair<coeffType, uint32_t> > > rightSide;
			
			/**
			 * The list of operations which have to be executed. One 'F4Operations' represents a set of operations which
			 * will be done parallel. No two elements of ops can be done in parallel, since there are dependencies or
			 * blocking operations.
			 */
			std::vector<F4Operations > ops;
			/**
			 * The dense representation of the matrix entries in the non pivot part. The operations of 'ops' will be executed
			 * on this matrix
			 */
			coeffMatrix matrix;

			/**
			 * For each row in the coefficient matrix this vector stores the information, how many rows have to be reduced
			 * before the given row is independent and can be used as operator.
			 */
			std::vector<size_t> deps;

			/**
			 * The number of rows in the matrix representing S-Polynomials.
			 */
			size_t upper;

			size_t reduceBlockSize;
			
			size_t rowCount;

			/**
			 * Stores during 'gauss' for each row if the given row is empty ( only zero entries )
			 */
			std::vector<bool> empty;

			/**
			 * Verbosity, which can be changed during runtime, nothing which should
			 * influence performance:
			 *
			 * 1 - Overall Runtime of the computation
			 * 2 - Sum of reduction time
			 * 4 - Sum of matrix setup/prepare time
			 * 8 - Sum of update time
			 * 16 - Print sugar degree during reduction step
			 * 32 - Print time of reduction step
			 * 64 - Print matrix size during reduction step 
			 *
			 * Add all debug informations to the verbosity flag of the ()-operator below, to get these
			 * details, e.g. to get all informations use '127' and to get only statistics at the end
			 * use '15'
			 */
			int verbosity;
			/**
			 * The output stream for the verbose messages. 
			 */
			std::ostream* out;

			/**
			 * Stores the time of all reduction steps
			 */
			double reductionTime;
			/**
			 * Stores the time of all matrix preparation steps
			 */
			double prepareTime;
			/**
			 * Stores the time of all critcial pairs update steps
			 */
			double updateTime;

			/**
			 * Parallel gaussian elemination on 'rs'. This function is called
			 * at the end of the reduction
			 */
			void gauss();

			/**
			 * Parallel reduction using all operations stored in 'ops'
			 */
			void pReduce();

			/**
			 * Setup the reduction matrix.
			 */
			void prepare();

			/**
			 * Do the next reduction step based on the given critical pairs.
			 * 'polys' will be used to store the result of the reduction step.
			 */
			void reduce(std::vector<Polynomial>& polys);

			/**
			 * Update the set of critical pairs using the polynomials 'polys',
			 * which are the new candidates for elements of the groebner basis.
			 */
			void updatePairs(std::vector<Polynomial>& polys, bool initial = false);
	
			/**
			 * Process a S-Polynomial (ir*current) and do the symbolic preprocessing by adding required
			 * reduction polynomials to the matrix
			 */
			void setupRow(Polynomial& current, Term& ir, size_t i, tbb::blocked_range<size_t>& range);

			/**
			 * Convert the non pivot part of the matrix from sparse ('rightSide') to dense ('rs') representation
			 */
			void setupDenseRow(coeffMatrix& rs, size_t offset, tbb::blocked_range<size_t>& range);
			
			/**
			 * Reduce a slice of the operation set ops[i] 
			 */
			void pReduceRange(coeffMatrix& rs, std::vector<size_t>& prefixes, std::vector<size_t>& suffixes, size_t i, tbb::blocked_range<size_t>& range);

			/**
			 * This operator computes a groebner basis for the given set of 'generators' using the term ordering 'O'.
			 * The polynomials should already be brought into the coefficient field 'field'. The number of threads
			 * is a upper bound for the parallelism and should be set to the number of available processor cores.
			 * 'verbosity' and 'output' can be used for debugging and performance measurement.
			 */
			std::vector<Polynomial> operator()(std::vector<Polynomial>& generators, const TOrdering* O, CoeffField* field, int threads = 1, int verbosity = 0, std::ostream& output = std::cout);
	};

	/**
	 * Helper class for parallel matrix setup. Will be used by tbb::parallel_for().
	 * The operator() is just a callback for the setupRow() function of the class F4
	 */
	struct F4SetupRow
	{
		/**
		 * The used instance of the F4 algorithm
		 */
		F4& f4;
		
		/**
		 * the polynomial which will be used to compute the S-Polynomial (current*ir)
		 */
		Polynomial& current;
		/**
		 * the term which will be used to compute the S-polynomial (current*ir)
		 */
		
		Term& ir;
		/**
		 * The index of the current row of the coefficient matrix
		 */
		size_t i;

		/**
		 * Construct a new instance of F4SetupRow.
		 */
		F4SetupRow(F4& f4, Polynomial& current, Term& ir, size_t i) : f4(f4), current(current), ir(ir), i(i) {}

		/**
		 * Call back the setupRow function of the given f4 instance.
		 */
		void operator() (tbb::blocked_range<size_t>& range) const
		{
			f4.setupRow(current, ir, i, range);
		}
	};

	/**
	 * Helper class for parallel matrix setup. Will be used by tb::parallel_for()
	 * The operator() is just a callback for the setupDenseRow() function of the class f4
	 */
	struct F4SetupDenseRow {
		/**
		 * The used instance of the F4 algorithm
		 */
		F4& f4;
		coeffMatrix& rs;
		size_t offset;

		/**
		 * Construct a new instance of F4SetupDenseRow
		 */
		F4SetupDenseRow(F4& f4, coeffMatrix& rs, size_t offset) : f4(f4), rs(rs), offset(offset) {}

		/**
		 * Call back the setupDenseRow function of the given f4 instance
		 */
		void operator() (tbb::blocked_range<size_t>& range) const
		{
			f4.setupDenseRow(rs, offset, range);
		}
	};

	/**
	 * Helper class for parallel matrix reduction. Will be used by tbb:parallel_for()
	 * The operator() is just a callback for the pReduceRange() function of the class f4
	 */
	struct F4PReduceRange
	{
		/**
		 * The used instance of the F4 algorithm
		 */
		F4& f4;
		coeffMatrix& rs;
		std::vector<size_t>& prefixes;
		std::vector<size_t>& suffixes;
	
		/**
		 * The index for the current set of operations (f4.ops[i])
		 */
		size_t i;

		/**
		 * Construct a new instance of F4PReduceRange
		 */
		F4PReduceRange(F4& f4, coeffMatrix& rs, std::vector<size_t>& prefixes, std::vector<size_t>& suffixes, size_t i) : f4(f4), rs(rs), prefixes(prefixes), suffixes(suffixes), i(i) {}

		/**
		 * Call back the pReduceRange function of the given f4 instance
		 */
		void operator() (tbb::blocked_range<size_t>& range) const
		{
			f4.pReduceRange(rs, prefixes, suffixes, i, range);
		}
	};
}
#endif
