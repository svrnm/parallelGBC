/**
 * Implementation of the CoeffField class, see include/CoeffField.H
 *
 ***
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
#include "../include/CoeffField.H"
#include <iostream>

CoeffField::CoeffField(coeffType modn) : modn(modn)
{
	// If SSE is enabled this vector is needed to do computations in the
	// coefficient field.

	// Preassign exps, logs and invs. The
	// following initalization code is inspired by Singular (kernel/modulop.cc, function npInitChar)
	exps.assign(modn, 0);
	logs.assign(modn, 0);
	invs.assign(modn, 0);

	logs[0] = 0;
	exps[0] = 1;

	// The following procedure only works if modn is greater than 2
	if(modn > 2) {
		coeffType i;
		coeffType w = 1;
		do
		{
			logs[1] = 0;
			w++;
			i = 0;
			do
			{
				i++;
				exps[i] = (w * exps[i-1]) % modn;
				logs[exps[i]] = i;
			} while( exps[i] != 1);
		} while( i != modn - 1);

		// Having a double length list of 'exps' allows to remove
		// the length check if multiplication is done, since by doing
		// exps[logs[a]+logs[b]] the sum of logs[a]+logs[b] is less
		// than two times the length of exps.
		exps.insert(exps.end(), exps.begin()+1, exps.end());

		// Precompute the inverse map
		for(coeffType i = 0; i < modn; i++) {
			invs[i] = exps[modn - 1 - logs[ i ] ];
		}
	} else {
		// Alternate setup for the case modn = 2
		logs[1] = 1;
		invs[1] = 1;
		exps.push_back( 1 );
	}
}

// Helper function for the setup of the SSE vectors, this
// means multiply o[k+d] with the logarithm of lc.
#define omulc(d) ( o[k+d] != 0 ? exps[o[k+d] + lc] : 0)

void CoeffField::mulSub(coeffRow& t, coeffRow& o, coeffType c, size_t prefix, size_t suffix) const {

// Check if SSE is enabled at compile time
#if PGBC_USE_SSE == 1
	// Precompute the logarithm of logs[c] to speed up the computation, since the lookup has only to be done once
	coeffType lc = logs[c];
	coeffRow oc(o.begin(), o.end());
	for(size_t i = prefix; i < suffix; i++) {
		oc[i] = o[i] != 0 ? exps[o[i] + lc] : 0;
	}

	// Set the vector x to the beginning of the target. This will load __COEFF_FIELD_INTVECSIZE values into the vector
	__m128i* x = (__m128i*) &(t[0]);
	__m128i* y = (__m128i*) &(oc[0]);
	__m128i modnvec = __COEFF_FIELD_VECSET1( modn );
	
	// set the first vector position by the given prefix (=0-padding at the beginning)
	prefix = prefix / __COEFF_FIELD_INTVECSIZE;
	suffix = suffix / __COEFF_FIELD_INTVECSIZE;
	// Iterate over the target+operator vector by doing __COEFF_FIELD_INTVECSIZE steps in parallel
	for(size_t i = prefix; i < suffix; i++ ) {
		// Read this as: x[i] = x[i] + ( (y > x[i] & modn) ) - y;
		//
		// This computes for __COEFF_FIELD_INTVECSIZE values the reduction at once. The reduction is done
		// by the following steps:
		// 1) y > x[1] = tmp0: Check if y is greater than x[i]. This returns a vector, which stores 0xf...f if an element of y is greater than x[i]
		// and 0x0...0 if not, where the length of the result matches the number of coefficient bits.
		// 2) (tmp0 & modn) = tmp1: For each position in y > x[i] step 1 has computed 0xf...f as result, so this operation
		// stores "modn" in this case in the result
		// 3) x[i] + tmp1 = tmp2: For each position in x[i] the value is increased by modn if the condition in 1 is true 
		// 4) x[i] + y = x[i]: Finally x[i] and y can be added.
		//
		x[i] = __COEFF_FIELD_VECADD(x[i], __COEFF_FIELD_VECSUB(__COEFF_FIELD_VECAND(__COEFF_FIELD_VECGT(y[i], x[i]), modnvec), y[i]));
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
			// Do the addition: If b is greater than t[k] an additional modn is
			// added to the result to stay within [0;modn[ at the end.
			t[k] = (b > t[k]) ? t[k] - b + modn : t[k] - b;
		}
	}
#endif
}
