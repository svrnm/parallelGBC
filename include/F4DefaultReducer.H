/**
 *  This file includes the headers for the 'defaultReducer' which implements the algorithm as described
 *  in http://link.springer.com/chapter/10.1007/978-3-642-32973-9_22 
 *
 ***********************
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
#ifndef F4_DEFAULTREDUCER_H
#define F4_DEFAULTREDUCER_H
#include <set>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/bind.hpp>
#include <tbb/concurrent_unordered_map.h>
#include <unordered_map>
#include <tbb/concurrent_vector.h>
#include <tbb/blocked_range.h>
#include <tbb/atomic.h>
#include "../include/Term.H"
#include "../include/Polynomial.H"
#include "../include/F4Reducer.H"
#include "../include/F4Algorithm.H"
#include "../include/F4Simplify.H"
#include "../include/F4SimplifyDB.H"

namespace parallelGBC {

	class F4Logger;

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


    template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {   
        ar & targets;
        ar & opers;
        ar & factors;
    }
	};

	class F4DefaultReducer : public virtual F4Reducer {
		public:
			int doSimplify;
			size_t reduceBlockSize;

			std::vector<std::pair<uint32_t, uint32_t> > newPivots;

			
			// this will become the 'real' part
			F4SimplifyDB* simplifyDB;
			std::vector<std::pair<size_t, Term> > rowOriginDB;
			
			// remains to make comparison possible
			F4Simplify* simplify;
			std::vector<std::pair<Term, Polynomial> > rowOrigin;

			std::vector<std::vector<std::pair<uint32_t, coeffType> > > savedRows; 
			std::vector<uint32_t> termMapping;

			/**
			 * The set of terms which occure during reduction and are not leading terms
			 */
			std::map<Term, uint32_t, Term::comparator> terms;

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

			tbb::concurrent_vector<std::pair<size_t, Term> > rows;

			/**
			 * For each row in the coefficient matrix this vector stores the information, how many rows have to be reduced
			 * before the given row is independent and can be used as operator.
			 */
			std::vector<size_t> deps;

			/**
			 * The number of rows in the matrix representing S-Polynomials.
			 */
			size_t upper;


			size_t rowCount;
			/**
			 * Stores during 'gauss' for each row if the given row is empty ( only zero entries )
			 */
			std::vector<bool> empty;

#if PGBC_WITH_MPI == 1
      std::vector<tbb::concurrent_vector<std::pair<uint32_t, std::pair<coeffType, uint32_t> > > > toSend;
			std::vector<std::vector<std::pair<uint32_t, std::pair<coeffType, uint32_t> > > > toSendCopy; // required :-(
#endif

			F4DefaultReducer(F4* f4, int doSimplify = 0, int reduceBlockSize = 1024) : F4Reducer(f4), doSimplify(doSimplify), reduceBlockSize(reduceBlockSize) {
				termCounter = 0;
				upper = 0;
				terms = std::map<Term, uint32_t, Term::comparator>( Term::comparator(f4->O, true) );
				pivotsOrdered = std::map<Term, uint32_t, Term::comparator>( Term::comparator(f4->O, true) );
				if(doSimplify == 2) {
					simplifyDB = new F4SimplifyDB(f4->O);
				} else if(doSimplify == 1) {
					simplify = new F4Simplify();
				}
			}

			~F4DefaultReducer() {
				if(doSimplify == 2) {
					delete simplifyDB;
				} else if(doSimplify == 1) {
					delete simplify;
				}
			}



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

			virtual void init() {
				if(doSimplify == 2) {
					for(size_t i = 0; i < f4->groebnerBasis.size(); i++) {
						Term one =  f4->groebnerBasis[i].LT().getOne();
						simplifyDB->insert(i, one, f4->groebnerBasis[i]);
					}
				}
			}

			virtual void reduce(std::vector<Polynomial>& polys, degreeType currentDegree);

			virtual void addSPolynomial(size_t i, size_t j, Term& lcm);

			virtual void client();

			virtual void finish();

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
			void pReduceRange(coeffMatrix& rs, std::vector<size_t>& prefixes, std::vector<size_t>& suffixes, size_t i, size_t offset, tbb::blocked_range<size_t>& range);

			void prepareOperator(coeffRow& row, size_t index, size_t& prefix, size_t& suffix, size_t offset); 
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
		F4DefaultReducer& reducer;

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
		F4SetupRow(F4DefaultReducer& reducer, Polynomial& current, Term& ir, size_t i) : reducer(reducer), current(current), ir(ir), i(i) {}

		/**
		 * Call back the setupRow function of the given f4 instance.
		 */
		void operator() (tbb::blocked_range<size_t>& range) const { reducer.setupRow(current, ir, i, range); }
	};

	/**
	 * Helper class for parallel matrix setup. Will be used by tb::parallel_for()
	 * The operator() is just a callback for the setupDenseRow() function of the class f4
	 */
	struct F4SetupDenseRow {
		/**
		 * The used instance of the F4 algorithm
		 */
		F4DefaultReducer& reducer;
		coeffMatrix& rs;
		size_t offset;

		/**
		 * Construct a new instance of F4SetupDenseRow
		 */
		F4SetupDenseRow(F4DefaultReducer& reducer, coeffMatrix& rs, size_t offset) : reducer(reducer), rs(rs), offset(offset) {}

		/**
		 * Call back the setupDenseRow function of the given f4 instance
		 */
		void operator() (tbb::blocked_range<size_t>& range) const { reducer.setupDenseRow(rs, offset, range); }
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
		F4DefaultReducer& reducer;
		coeffMatrix& rs;
		std::vector<size_t>& prefixes;
		std::vector<size_t>& suffixes;

		/**
		 * The index for the current set of operations (f4.ops[i])
		 */
		size_t i;
		size_t offset;

		/**
		 * Construct a new instance of F4PReduceRange
		 */
		F4PReduceRange(F4DefaultReducer& reducer, coeffMatrix& rs, std::vector<size_t>& prefixes, std::vector<size_t>& suffixes, size_t i, size_t offset) : reducer(reducer), rs(rs), prefixes(prefixes), suffixes(suffixes), i(i), offset(offset) {}

		/**
		 * Call back the pReduceRange function of the given reducer instance
		 */
		void operator() (tbb::blocked_range<size_t>& range) const { reducer.pReduceRange(rs, prefixes, suffixes, i, offset, range); }
	};
}
#endif
