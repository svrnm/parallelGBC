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
#ifndef T_MONOID_H
#define T_MONOID_H
#include <stdint.h>
#include <string>
#include <vector>
#include <unordered_set>
#include "F4Utils.H"

using namespace std;

#ifndef F4_DEGREE_BITS
#ifdef F4_COEFF_BITS
#define F4_DEGREE_BITS F4_COEFF_BITS
#else
#define F4_DEGREE_BITS 32
#endif
#endif

#if F4_DEGREE_BITS <= 16
#define __MONOID_DEGREE_TYPE int16_t
#else
#define __MONOID_DEGREE_TYPE int32_t
#endif


typedef __MONOID_DEGREE_TYPE degreeType;

class Term;

class TMonoid {
	protected:
		struct TermEquals : std::binary_function<const Term* const, const Term* const , bool>
	{
		bool operator()(const Term* const t1, const Term* const t2) const ;
	};
		struct TermHash : std::unary_function<const Term* const, std::size_t>
	{
		std::size_t operator()(const Term* const t ) const;
	};

		typedef unordered_set<Term*, TermHash, TermEquals> TermSet;
		TermSet terms;
		Term* one;

	public:
		const size_t N;
		const size_t D;
		TMonoid(size_t N);
		

		const Term* createElement(const string& s, degreeType min);
		const Term* createElement(vector<long>& v);
		const Term* createElement(Term* t);

		size_t size() const {
			return terms.size();
		}
		const Term* getOne() const {
			return this->one;
		}
		~TMonoid();
	private:
		TMonoid(const TMonoid& m) : N(m.N), D(m.D) { }
};
#endif
