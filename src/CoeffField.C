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

#ifdef __SSE2__
void CoeffField::mulSub(std::vector<coeffType>& t, std::vector<coeffType>& o, coeffType c, size_t prefix, size_t suffix) const {
	std::vector<coeffType> tmp(o);
	if(c != 1) {
		coeffType lc = logs[c];
		for(size_t k = prefix; k < suffix; k++) {
			if(tmp[k] != 0) {
				tmp[k] = exps[logs[tmp[k]] + lc];
			}
		}
	}

	__m128i* y = (__m128i*) &(tmp[0]);
	__m128i* x = (__m128i*) &(t[0]);
	size_t n = suffix / (__COEFF_FIELD_INTVECSIZE);
	for(size_t k = prefix / __COEFF_FIELD_INTVECSIZE; k < n; k++) {
		// x[k] + = ( (y[k] > x[k]) & modn) - y[k];
		x[k] = __COEFF_FIELD_VECADD(x[k], __COEFF_FIELD_VECSUB(__COEFF_FIELD_VECAND(__COEFF_FIELD_VECGT(y[k], x[k]), modnvec), y[k]));
	}
}
#else
void CoeffField::mulSub(std::vector<coeffType>& t, std::vector<coeffType>& o, coeffType c, size_t prefix, size_t suffix) const {
	c = logs[c];
	for(size_t k = prefix; k < suffix; k++) {
		if(o[k] != 0) {
			coeffType b = exps[logs[o[k]] + c];
			t[k] = (b > t[k]) ? t[k] - b + modn : t[k] - b;
		}
	}
}
#endif
