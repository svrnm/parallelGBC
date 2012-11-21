/**
 * This file includes the headers of 'TMonoid' which represents the monoid containing all instances
 * of terms with a given number of indetermenants. The monoid is just a container object for the
 * instances of terms. The monoid has a hash map, which stores for each Term a unique instance. 
 * This makes it possible that a Term is only identified by a small amount of memory (the pointer
 * to the instance) and comparision of equality is just a comparision of pointers. 
 *
 * The details of TermInstance* and Term are given in the header file include/Term.h
 *
 *******
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
#ifndef T_MONOID_H
#define T_MONOID_H
#include <stdint.h>
#include <string>
#include <vector>
#include <unordered_set>
#include <tbb/concurrent_unordered_set.h>

/**
 * The number of 'bits' can be defined independently of
 * the number of 'bits' for the coefficients. Normally
 * these two values should be equal, therefor the default
 * behaviour is to set the degree bits equal to the coefficient
 * bits
 */
#ifndef PGBC_DEGREE_BITS
	#ifdef PGBC_COEFF_BITS
		#define PGBC_DEGREE_BITS PGBC_COEFF_BITS
	#else
		#define PGBC_DEGREE_BITS 32
	#endif
#endif

/**
 * Choose the degreeType, you can choose between 16 and 32bit
 * TODO: Why not 8 bit?
 */
#if PGBC_DEGREE_BITS <= 16
	typedef int16_t degreeType;
#else
	typedef int32_t degreeType;
#endif

/**
 * Forward declaration of TermInstance. This class is realized in include/Term.H
 */
class TermInstance;
/**
 * Forward declaration of TOrdering. This class is realized in include/TOrderin.H
 */
class TOrdering;

/**
 * The class TMonoid represents a "Term Monoid" which is the set of all terms of
 * degree 0 to N. The given implementation of the monoid is only a storage, the
 * computation is mostly done directly by the class TermInstance.
 *
 * The TMonoid is a hash map which takes care that there are within the monoid
 * unique instances of terms, which makes it possible to save a term as a single
 * pointer to this object.
 */
class TMonoid {
	protected:
		/**
		 * Helper class, to compare two instances of TermInstance for equality
		 */
		struct TermInstanceEquals : std::binary_function<const TermInstance* const, const TermInstance* const , bool>
		{
			bool operator()(const TermInstance* const t1, const TermInstance* const t2) const ;
		};
		/**
		 * Helper class, to retrieve the nash of a given TermInstance.
		 */
		struct TermInstanceHash : std::unary_function<const TermInstance* const, std::size_t>
		{
			std::size_t operator()(const TermInstance* const t ) const;
		};

		/**
		 * Shortname for the set of all term instances container type
		 */
		typedef tbb::concurrent_unordered_set<TermInstance*, TermInstanceHash, TermInstanceEquals> TermInstanceSet;

		/**
		 * The set which contains all unique copies of TermInstances
		 */
		TermInstanceSet terms;
		/**
		 * An excplicit pointer to the element "1", i.e. the term with degree 0
		 */
		TermInstance* one;

	public:
		/**
		 * The number of indeterminants
		 */
		const size_t N;
		
		/**
		 * A default copy of the lexicographic term ordering
		 */
		TOrdering* lexOrdering;
		/**
		 * A default copy of the degree lexicographic term ordering
		 */
		TOrdering* degLexOrdering;

		/**
		 * Construct a new instance of a TMonoid with N indeterminants
		 */
		TMonoid(size_t N);

		/**
		 * Return a pointer to the unique 'TermInstance' which is equivalent to the
		 * given representation by the string 's'. The parameter 'min' sets the smallest
		 * indeterminant. This allows to read in x[1] as the smallest indeterminant, while
		 * internal the index 0 is used.
		 *
		 * The format of 's':
		 * - An indeterminant is given by x[i]
		 * - The power of an indeterminant is given by x[i]^n
		 * - The product of several powers of indeterminants is given by x[i]^n*x[j]^m*...
		 */
		const TermInstance* createElement(const std::string& s, degreeType min);
		
		/**
		 * Return a pointer to the unique 'TermInstance' which is equivalent to the
		 * given representation by the degree vector 'v'
		 */
		const TermInstance* createElement(const std::vector<degreeType>& v);

		/**
		 * Check if an other 'TermInstance' which is equivalent to t already exists. If
		 * so, t is deleted and the pointer to the existing instance is given. If not,
		 * t is inserted into the monoid and returned.
		 */
		const TermInstance* createElement(TermInstance* t);

		/**
		 * Return a pointer to the lexicographic ordering
		 */
		TOrdering* lex();

		/**
		 * Return a pointer to the degree lexicographic ordering
		 */
		TOrdering* degLex();
		
		/**
		 * Returns the number of terms, which are in use. This is not
		 * the size of the mathematical object.
		 */
		size_t size() const {
			return terms.size();
		}

		/**
		 * Returns a pointer to the term '1'
		 */
		const TermInstance* getOne() const {
			return this->one;
		}

		/**
		 * Destruct the TMonoid object.
		 */
		~TMonoid();

	private:
		/**
		 * Private copy constructor, should never be used, since a copy
		 * would destroy the uniqueness of the already existing terms.
		 */
		TMonoid(const TMonoid& m) : N(m.N) { }
};
#endif
