/*
 *  Header file for the coefficient field. As coefficient field you can use any
 *  prime field within 2<=p<=2^k, where p is prime and k = PGBC_COEFF_BITS. The
 *  default size of PGBC_COEFF_BITS is 16, so you can use e.g. the field with 32003
 *  elements.
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
#ifndef T_COEFF_FIELD_H
#define T_COEFF_FIELD_H
#include <stdint.h>
#include <vector>
#include <iostream>
#include <xmmintrin.h>
#include <mmintrin.h>
#define BREAKPOINT {while(getchar() != '\n');}

// By default we use 16 bits for the coefficients
#ifndef PGBC_COEFF_BITS 
	#define PGBC_COEFF_BITS 16
#endif

// Define the SSE functions in dependency of PGBC_COEFF_BITS
#if PGBC_COEFF_BITS <= 8
	#define __COEFF_FIELD_INTVECSIZE 16
	#define __COEFF_FIELD_VECGT(x,y) (_mm_cmpgt_epi8(x,y))
	#define __COEFF_FIELD_VECSET1(x) (_mm_set1_epi8(x))
	#define __COEFF_FIELD_VECADD(x,y) (_mm_add_epi8(x,y))
	#define __COEFF_FIELD_VECSUB(x,y) (_mm_sub_epi8(x,y))
	#define __COEFF_FIELD_NUMBERTYPE int8_t
#else
	#if PGBC_COEFF_BITS <= 16
		#define __COEFF_FIELD_INTVECSIZE 8
		#define __COEFF_FIELD_VECGT(x,y) (_mm_cmpgt_epi16(x,y))
		#define __COEFF_FIELD_VECSET1(x) (_mm_set1_epi16(x))
		#define __COEFF_FIELD_VECADD(x,y) (_mm_add_epi16(x,y))
		#define __COEFF_FIELD_VECSUB(x,y) (_mm_sub_epi16(x,y))
		#define __COEFF_FIELD_NUMBERTYPE int16_t
	#else
		#define __COEFF_FIELD_INTVECSIZE 4
		#define __COEFF_FIELD_VECGT(x,y) (_mm_cmpgt_epi32(x,y))
		#define __COEFF_FIELD_VECSET1(x) (_mm_set1_epi32(x))
		#define __COEFF_FIELD_VECADD(x,y) (_mm_add_epi32(x,y))
		#define __COEFF_FIELD_VECSUB(x,y) (_mm_sub_epi32(x,y))
		#define __COEFF_FIELD_NUMBERTYPE int32_t
	#endif
#endif

#define __COEFF_FIELD_VECAND(x,y) (_mm_and_si128(x,y))

// define the type for coefficients
typedef __COEFF_FIELD_NUMBERTYPE coeffType;

/*
 * This class represents the coefficient field. A coefficient is a
 * primitive type (i.e. coeffType = intXX_t) and therefore the functions
 * are applied via an instance of CoeffField. This is good for performance
 * and memory usage, not for software design ...
 */
class CoeffField
{
	public:
		// The modulus, i.e. n of F_n
		static const size_t pad = __COEFF_FIELD_INTVECSIZE * 4;
		const coeffType modn;

	private:
		// precalculated exponents and logarithms
		std::vector<coeffType> exps;
		std::vector<coeffType> logs;
		std::vector<coeffType> invs;

		// array with k times the modules for SSE computations
		__m128i modnvec;

		// this property is used by mulSub as multiplication factor
		coeffType factor;

		// two helpers for the SSE computations in mulSub
		size_t limit;
		size_t s;

	public:
		CoeffField(coeffType modn);
		coeffType mul(coeffType a, coeffType b) const
		{ 
			if( a == 0 || b == 0 ) return 0;
			return exps[logs[a] + logs[b]];
		}

		coeffType bringIn(coeffType a) const
		{
			return ((a % modn) + modn) % modn;
		}

		coeffType inv(coeffType a) const
		{ 
			return invs[a];
		}
		coeffType minus(coeffType a) const
		{
			if( a == 0 ) return 0;
			return modn - a;
		}
		coeffType sub(coeffType a, coeffType b) const
		{
			if( b == 0) return a;
			return b > a ? a + modn - b: a - b;
		}
		coeffType add(coeffType a, coeffType b) const
		{
			coeffType c = a + b;
			return c < modn ? c : c - modn;
		}

		coeffType getFactor(coeffType a) {
			return logs[a];
		}

		coeffType mulSub(coeffType a, coeffType b, coeffType c) const {
			if(b == 0) { return a; }
			b = exps[logs[b] + c];
			return b > a ? a + modn - b : a - b;
		}

		void mulSub(std::vector<coeffType>& t, std::vector<coeffType>& o, coeffType c, size_t prefix, size_t suffix) const;
};
#endif
