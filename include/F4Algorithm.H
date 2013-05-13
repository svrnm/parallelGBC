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
#include <tbb/concurrent_unordered_map.h>
#include <unordered_map>
#include <tbb/concurrent_vector.h>
#include <tbb/blocked_range.h>
#include <tbb/atomic.h>
#include "../include/Term.H"
#include "../include/Polynomial.H"
#include "../include/F4Utils.H"
#include "../include/F4Logger.H"
#include "../include/F4Reducer.H"
#if PGBC_WITH_MPI == 1
#include <boost/mpi.hpp>
#include <boost/serialization/utility.hpp>
#endif

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

			bool withSugar;

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
			
			F4Logger* log;
		
			F4Reducer* reducer;

			/**
			 * Update the set of critical pairs using the polynomials 'polys',
			 * which are the new candidates for elements of the groebner basis.
			 */
			void updatePairs(std::vector<Polynomial>& polys, bool initial = false);

			void select();
			
			void reduce(std::vector<Polynomial>& polys);

			void setReducer(F4Reducer* reducer) {
				this->reducer = reducer;
			}

#if PGBC_WITH_MPI == 1

			boost::mpi::communicator& world;
			tbb::concurrent_vector<boost::mpi::request> reqs;

			F4(const TOrdering* O, CoeffField* field, boost::mpi::communicator& world, bool withSugar = true, int threads = 1, int verbosity = 0, std::ostream& output = std::cout) : O(O), field(field), withSugar(withSugar), threads(threads), log(new F4Logger(verbosity, &output)), world(world) { }

#else 
			F4(const TOrdering* O, CoeffField* field, bool withSugar = true, int threads = 1, int verbosity = 0, std::ostream& output = std::cout) : O(O), field(field), withSugar(withSugar), threads(threads), log(new F4Logger(verbosity, &output)) { }
#endif


			~F4() {
				delete log;
			}

			/**
			 * This operator computes a groebner basis for the given set of 'generators' using the term ordering 'O'.
			 * The polynomials should already be brought into the coefficient field 'field'. The number of threads
			 * is a upper bound for the parallelism and should be set to the number of available processor cores.
			 * 'verbosity' and 'output' can be used for debugging and performance measurement.
			 */
			std::vector<Polynomial> compute(std::vector<Polynomial>& generators);
	};
}
#endif
