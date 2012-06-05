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

class TermInstance;

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

		typedef unordered_set<TermInstance*, TermInstanceHash, TermInstanceEquals> TermInstanceSet;
		TermInstanceSet terms;
		TermInstance* one;

	public:
		const size_t N;
		const size_t D;
		TMonoid(size_t N);
		

		const TermInstance* createElement(const string& s, degreeType min);
		const TermInstance* createElement(const vector<degreeType>& v);
		const TermInstance* createElement(TermInstance* t);

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
