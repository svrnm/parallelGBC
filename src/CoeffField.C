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
#include "../include/CoeffField.H"
#include <iostream>

CoeffField::CoeffField(coeffType modn) : modn(modn)
{
	modnvec = __COEFF_FIELD_VECSET1( modn );

	coeffType i;
	coeffType w = 1;
	exps.assign(modn, 0);
	logs.assign(modn, 0);
	invs.assign(modn, 0);
	logs[0] = 0;
	exps[0] = 1;
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
	exps.insert(exps.end(), exps.begin()+1, exps.end());


	for(coeffType i = 0; i < modn; i++) {
		invs[i] = exps[modn - 1 - logs[ i ] ];
	}
}

#define omulc(d) (o[k+d] != 0 ? exps[o[k+d] + lc] : 0)

#if PGBC_USE_SSE == 1
void CoeffField::mulSub(std::vector<coeffType>& t, std::vector<coeffType>& o, coeffType c, size_t prefix, size_t suffix) const {
	__m128i* x = (__m128i*) &(t[0]);
	coeffType lc = logs[c];
	size_t i = prefix / __COEFF_FIELD_INTVECSIZE;
	for(size_t k = prefix; k < suffix; k+=__COEFF_FIELD_INTVECSIZE) {
		// Not a beauty ...
		#if PGBC_COEFF_BITS <= 8
			__m128i y = _mm_set_epi8(omulc(15),omulc(14),omulc(13),omulc(12),omulc(11),omulc(10),omulc(9),omulc(8), omulc(7),omulc(6),omulc(5),omulc(4),omulc(3),omulc(2),omulc(1),omulc(0));
		#else
			#if PGBC_COEFF_BITS <= 16
				__m128i y = _mm_set_epi16(omulc(7),omulc(6),omulc(5),omulc(4),omulc(3),omulc(2),omulc(1),omulc(0));
			#else
				__m128i y = _mm_set_epi32(omulc(3),omulc(2),omulc(1),omulc(0));
			#endif
		#endif
		// x[i] = x[i] + ( (y > x[i] & modn) ) - y;
		x[i] = __COEFF_FIELD_VECADD(x[i], __COEFF_FIELD_VECSUB(__COEFF_FIELD_VECAND(__COEFF_FIELD_VECGT(y, x[i]), modnvec), y));
		i++;
	}
}
#else
void CoeffField::mulSub(std::vector<coeffType>& t, std::vector<coeffType>& o, coeffType c, size_t prefix, size_t suffix) const {
	c = logs[c];
	for(size_t k = prefix; k < suffix; k++) {
		if(o[k] != 0) {
			coeffType b = exps[o[k] + c];
			t[k] = (b > t[k]) ? t[k] - b + modn : t[k] - b;
		}
	}
}
#endif
