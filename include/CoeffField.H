/**
 *  Header file for the coefficient field. As coefficient field you can use any
 *  prime field within 2<=p<=2^k, where p is prime and k = PGBC_COEFF_BITS. The
 *  default size of PGBC_COEFF_BITS is 16, so you can use e.g. the field with 32003
 *  elements.
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
#ifndef T_COEFF_FIELD_H
#define T_COEFF_FIELD_H
#include <stdint.h>
#include <vector>
#include <iostream>
#include <xmmintrin.h>
#include <mmintrin.h>
#include <tbb/concurrent_vector.h>

/** 
 * Coefficents in the finite field are represented by 
 * primitive data types (like int16_t), so there is no
 * class like Coefficient. The constand PGBC_COEFF_BITS
 * can be set at compile time to define the bitsize of
 * the primitive types, default is 16.
 */
#ifndef PGBC_COEFF_BITS 
#define PGBC_COEFF_BITS 16
#endif

/**
 * The SSE implementation of mulSub(coeffRow...) requires a certain set of functions
 * which are generalized here since the called functions are dependend of the
 * used primitive data type.
 *
 * The INTVECSIZE represents the vector size and the larger this value is, the
 * more operations can be done simultaneously. The NUMBERTYPE abstracts the
 * primitive data type. The number of available bits defines a upper bound
 * for the modulus of the finite field.
 *
 * The possible primitive data types are:
 *
 * int8_t: 8-bit per coefficient, so vector size is 16
 * int16_t: 16-bit per coefficient, so vector size is 8
 * int32_t: 32-bit per coefficient, so vector size is 4
 *
 * The following functions are required:
 * VECGT(x,y): Greater than (>) for all elements in the sse vector
 * VECSET1(x): Setup a vector which contains n times the element x
 * VECADD(x,y): Vector addition (+)
 * VECSUB(x,y): Vector subtraction (-)
 *
 * PGBC_COEFF_BITS is defining the size. The user should use a value like 8,16 or
 * 32, all other types are rounded up to the next bigger value.
 */
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

/**
 * VECAND: Logical AND (&) for two 128-bit values. This function
 * is independent of the used NUMBERTYPE
 */
#define __COEFF_FIELD_VECAND(x,y) (_mm_and_si128(x,y))

/**
 * coeffType is the short name for the used coefficient type
 * coeffRow represents a vector of coefficients
 * coeffMatrix represents a vector of vectors of coefficients, therfore
 * a matrix of coefficients. 
 */
typedef __COEFF_FIELD_NUMBERTYPE coeffType;
typedef std::vector<coeffType> coeffRow;
typedef std::vector<coeffRow> coeffMatrix;

/*
 * This class represents the coefficient field. A coefficient is a
 * primitive type as defined above (i.e. coeffType = intXX_t). The functions
 * are applied via an instance of CoeffField.
 *
 * For all operations except bringIn you have to take care, that the operands
 * are smaller than the given modulus "modn".
 */
class CoeffField
{
	public:
		/**
		 * Padding for the vectors used in SSE mulSub(coeffRow...). This is required to make
		 * SSE operations available for vectors not having a matching size (i.e.
		 * divisible by __COEFF_FIELD_INTVECSIZE={4,8,16}
		 */
		static const size_t pad = __COEFF_FIELD_INTVECSIZE;

		/**
		 * The modulus, i.e. n of the finite field F_n
		 */
		const coeffType modn;

	protected:
		/**
		 * precalculated tables for exponents, logarithms and inverse
		 * Attention! All this vectors will have size modn or 2*modn, so
		 * if the NUMBERTYPE is for example int32_t, there can be consumed
		 * a lot of memory if modn is big.
		 */
		coeffRow exps;
		coeffRow logs;
		coeffRow invs;

		/**
		 * The SSE implementation of mulSub(coeffRow...) requires a vector which contains
		 * __COEFF_FIELD_INTVECSIZE times the value of modn.
		 */
		__m128i modnvec;

	public:
		/**
		 * Constructs and setups the coefficient field. 
		 * The parameter modn defines the size of the field, meaning that
		 * the field has the elements {0,...,modn-1}
		 *
		 * This method sets up the tables exps, logs and invs.
		 */
		CoeffField(coeffType modn);

		/**
		 * Multiplication of two coefficient types. The operation is done
		 * by doing exp(log(a) + log(b)), which is faster than a normal
		 * multiplication (*) and the modulus operation (%)
		 */
		coeffType mul(coeffType a, coeffType b) const
		{ 
			if( a == 0 || b == 0 ) return 0;
			return exps[logs[a] + logs[b]];
		}

		/**
		 * Division of two coefficient types. Assert that b is not equal zero,
		 * behaviour is not certain.
		 */
		coeffType div(coeffType a, coeffType b) const
		{
			return mul(a, inv(b));
		}

		/**
		 * Find a representation of a in the given finite field:
		 *
		 * (a % modn) may be negative if a is a negative.
		 * (a % modn) is smaller than modn, so (a % modn) + modn
		 * is a positive number. 
		 * If a was a positive number, (a % modn) + modn may be
		 * greater than modn, so ((a % modn) + modn) % modn is less
		 * than modn.
		 */
		coeffType bringIn(coeffType a) const
		{

			return ((a % modn) + modn) % modn;
		}

		/**
		 * Return the multiplicative inverse of a
		 */
		coeffType inv(coeffType a) const
		{ 
			return invs[a];
		}

		/**
		 * Return the additive inverse of a 
		 */
		coeffType minus(coeffType a) const
		{
			if( a == 0 ) return 0;
			return modn - a;
		}

		/**
		 * Subtract a from b and return the result.
		 *
		 * If b is zero, a-b = a.
		 * Else if b is greater than a the result would
		 * be negative. So return a + minus(b), else a-b.
		 */
		coeffType sub(coeffType a, coeffType b) const
		{
			if( b == 0) return a;
			return b > a ? a + modn - b: a - b;
		}

		/**
		 * Add a to b and return the result.
		 * if a+b is greater than modn, bring in c.
		 */
		coeffType add(coeffType a, coeffType b) const
		{
			coeffType c = a + b;
			return c < modn ? c : c - modn;
		}

		/**
		 * Return the logarithm of a. The function is called
		 * get factor, since it is used by the mulSub(...) functions
		 * for the third parameter "c" which is the factor of the
		 * matrix row operation "row a" minus "row b" * c.
		 */
		coeffType getFactor(coeffType a) {
			return logs[a];
		}

		/**
		 * Return a-b*exps[c]:
		 *
		 * This operation is used for elementwise reduction of matrix rows.
		 * The coefficient c is given as logarithm of the factor, since by
		 * doing this the table lookup has only to be done once, which
		 * improves performance measurable.
		 *
		 * TODO: This function is not used any more in the F4 implementation,
		 * maybe c can be changed again.
		 */
		coeffType mulSub(coeffType a, coeffType b, coeffType c) const {
			if(b == 0) { return a; }
			b = exps[logs[b] + c];
			return b > a ? a + modn - b : a - b;
		}


// Helper function for the setup of the SSE vectors, this
// means multiply o[k+d] with the logarithm of lc.
#define omulc(d) ( o[k+d] != 0 ? exps[o[k+d] + lc] : 0)



		/**
		 * Retunr a-b*exps[c] for all a in vector t and all b in vector o:
		 *
		 * This operation is used for the reduction of matrix row t. The additional parameter
		 * prefix and suffix are given for the row o, which might have a longer padding of 0s
		 * in the end and the beginning.
		 *
		 * The implementation is given in src/CoeffField.C and might be different if SSE is 
		 * enabled at compile time
		 */
		void mulSub(coeffRow& t, coeffRow& o, coeffType c, size_t prefix, size_t suffix) const {

			// Check if SSE is enabled at compile time
#if PGBC_USE_SSE == 1
			// Set the vector x to the beginning of the target. This will load __COEFF_FIELD_INTVECSIZE values into the vector
			__m128i* x = (__m128i*) &(t[0]);
			// Precompute the logarithm of logs[c] to speed up the computation, since the lookup has only to be done once
			coeffType lc = logs[c];
			// set the first vector position by the given prefix (=0-padding at the beginning)
			size_t i = prefix / __COEFF_FIELD_INTVECSIZE;
			// Iterate over the target+operator vector by doing __COEFF_FIELD_INTVECSIZE steps in parallel
			for(size_t k = prefix; k < suffix; k+=__COEFF_FIELD_INTVECSIZE) {
				// BEGIN: Setup the operator vector
				// LOAD __COEFF_FIELD_INTVECSIZE values into the operator y, the values of the vector o are multiplied with lc before
				// they are stored in the operator vector.
#if PGBC_COEFF_BITS <= 8
				__m128i y = _mm_set_epi8(omulc(15),omulc(14),omulc(13),omulc(12),omulc(11),omulc(10),omulc(9),omulc(8), omulc(7),omulc(6),omulc(5),omulc(4),omulc(3),omulc(2),omulc(1),omulc(0));
#else
#if PGBC_COEFF_BITS <= 16
				__m128i y = _mm_set_epi16(omulc(7),omulc(6),omulc(5),omulc(4),omulc(3),omulc(2),omulc(1),omulc(0));
#else
				__m128i y = _mm_set_epi32(omulc(3),omulc(2),omulc(1),omulc(0));
#endif
#endif
				// END: Setup the operator vector
				// Read this as: x[i] = x[i] + ( (y > x[i] & modn) ) - y;
				// This computes for __COEFF_FIELD_INTVECSIZE values the reduction at once. The reduction is done
				// by the following steps:
				// 1) y > x[1] = tmp0: Check if y is greater than x[i]. This returns a vector, which stores 0xf...f if an element of y is greater than x[i]
				// and 0x0...0 if not, where the length of the result matches the number of coefficient bits.
				// 2) (tmp0 & modn) = tmp1: For each position in y > x[i] step 1 has computed 0xf...f as result, so this operation
				// stores "modn" in this case in the result
				// 3) x[i] + tmp1 = tmp2: For each position in x[i] the value is increased by modn if the condition in 1 is true 
				// 4) x[i] + y = x[i]: Finally x[i] and y can be added.
				//
				x[i] = __COEFF_FIELD_VECADD(x[i], __COEFF_FIELD_VECSUB(__COEFF_FIELD_VECAND(__COEFF_FIELD_VECGT(y, x[i]), modnvec), y));
				i++;
			}
			// If SSE is not enabled at compile time
#else
			// Precompute the logarithm of logs[c] to speed up the computation, since the lookup has only to be done once
			c = logs[c];
			// Iterate from prefix to suffix and reduce the given elements in the target vector.
			// prefix and suffix are the zero-paddings of the operator vector.
			for(size_t k = prefix; k < suffix; k++) {
				// If o[k] is zero the result is t[k]
				if(o[k] != 0) {
					// Multiply o[k] with c
					coeffType b = exps[o[k] + c];
					// Do the addition: If b is greater than t[k] an additional modn is added to the result to stay within [0;modn[ at the end.
					t[k] = (b > t[k]) ? t[k] - b + modn : t[k] - b;
				}
			}
#endif
		}

};
#endif
