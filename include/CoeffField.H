/*
 *  Header file for the coefficient field. As coefficient field you can use any
 *  prime field within 2<=p<=2^k, where p is prime and k = F4_COEFF_BITS. The
 *  default size of F4_COEFF_BITS is 16, so you can use e.g. the field with 32003
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
#include <assert.h>
#include <stdint.h>
#include <vector>
#include <iostream>

// By default we use sse.
#ifndef F4_USE_SSE
	#define F4_USE_SSE 1
#endif

// By default we use 16 bits for the coefficients
#ifndef F4_COEFF_BITS 
	#define F4_COEFF_BITS 16
#endif

// Define some variables and functions in dependency of F4_COEFF_BITS
#if F4_COEFF_BITS <= 8
	#define __COEFF_FIELD_INTVECSIZE 16
	#define __COEFF_FIELD_GTVEC(x,y) (__builtin_ia32_pcmpgtb128 (x,y))
	#define __COEFF_FIELD_NUMBERTYPE char
#else
	#if F4_COEFF_BITS <= 16
		// 8*16 = 128 bits for the SSE vector.
		#define __COEFF_FIELD_INTVECSIZE 8
		// The build in GCC function to compare 8 16-bit integers at once
		#define __COEFF_FIELD_GTVEC(x,y) (__builtin_ia32_pcmpgtw128 (x,y))
		// The used type
		#define __COEFF_FIELD_NUMBERTYPE int16_t
	#else
		// 4*32 = 128 bits for the SSE vector
		#define __COEFF_FIELD_INTVECSIZE 4
		// The build in GCC function to compare 4 32-bit integers at once
		#define __COEFF_FIELD_GTVEC(x,y) (__builtin_ia32_pcmpgtd128 (x,y))
		// The used type
		#define __COEFF_FIELD_NUMBERTYPE int32_t
	#endif
#endif

using namespace std;

// define the type for coefficients
typedef __COEFF_FIELD_NUMBERTYPE coeffType;

// Define SSE Vector
typedef coeffType intvec __attribute__ ((vector_size (16)));

// Helper to initialize the vector
union intvecarray { intvec v; coeffType f[__COEFF_FIELD_INTVECSIZE]; };

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
		const coeffType modn;

	private:
		// precalculated exponents and logarithms
		vector<coeffType> exps;
		vector<coeffType> logs;
		vector<coeffType> invs;

		// modules minus 1
		const coeffType modnM1;

		// array with k times the modules for SSE computations
		intvecarray modnvec;

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

		#if F4_USE_SSE
		void mulSub(vector<coeffType>& t, vector<coeffType>& o, coeffType c) const {
			vector<coeffType> tmp(o);
			if(c != 1) {
				coeffType lc = logs[c];
				for(size_t k = 0; k < o.size(); k++) {
					if(o[k] != 0) { 
						if(o[k] == 1) {
							tmp[k] = c;
						} else {
							tmp[k] = exps[logs[o[k]] + lc]; 
						}
					}
				}
			}
			intvec* y = (intvec*) &(tmp[0]);
			intvec* x = (intvec*) &(t[0]);
		
			for(size_t k = 0; k < t.size() / (__COEFF_FIELD_INTVECSIZE); k++) {
				x[k+0] += (__COEFF_FIELD_GTVEC(y[k+0], x[k+0]) & modnvec.v) - y[k+0];
			}
			
			tmp.clear();
		}
		#else
		void mulSub(vector<coeffType>& t, vector<coeffType>& o, coeffType c) const {
			c = logs[c];
			for(size_t k = 0; k < o.size(); k++) {
				if(o[k] != 0) {
					coeffType b = exps[logs[o[k]] + c];
					t[k] = (b > t[k]) ? t[k] - b + modn : t[k] - b;	
				}
			}
		}
		#endif
};
#endif
