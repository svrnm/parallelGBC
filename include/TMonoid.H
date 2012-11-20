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

#ifndef PGBC_DEGREE_BITS
	#ifdef PGBC_COEFF_BITS
		#define PGBC_DEGREE_BITS PGBC_COEFF_BITS
	#else
		#define PGBC_DEGREE_BITS 32
	#endif
#endif

#if PGBC_DEGREE_BITS <= 16
	typedef int16_t degreeType;
#else
	typedef int32_t degreeType;
#endif

class TermInstance;
class TOrdering;

class TMonoid {
	protected:
		struct TermInstanceEquals : std::binary_function<const TermInstance* const, const TermInstance* const , bool>
		{
			bool operator()(const TermInstance* const t1, const TermInstance* const t2) const ;
		};
		struct TermInstanceHash : std::unary_function<const TermInstance* const, std::size_t>
		{
			std::size_t operator()(const TermInstance* const t ) const;
		};

		typedef tbb::concurrent_unordered_set<TermInstance*, TermInstanceHash, TermInstanceEquals> TermInstanceSet;

		TermInstanceSet terms;
		TermInstance* one;

	public:
		const size_t N;
		const size_t D;
		TOrdering* lexOrdering;
		TOrdering* degLexOrdering;
		TMonoid(size_t N);

		const TermInstance* createElement(const std::string& s, degreeType min);
		const TermInstance* createElement(const std::vector<degreeType>& v);
		const TermInstance* createElement(TermInstance* t);
		TOrdering* lex();
		TOrdering* degLex();
		

		size_t size() const {
			return terms.size();
		}
		const TermInstance* getOne() const {
			return this->one;
		}
		~TMonoid();
	private:
		TMonoid(const TMonoid& m) : N(m.N), D(m.D) { }
};
#endif
