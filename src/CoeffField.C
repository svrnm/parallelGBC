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

CoeffField::CoeffField(coeffType modn) : modn(modn), modnM1(modn-1)
{
	for(size_t k = 0; k < __COEFF_FIELD_INTVECSIZE; k++) {
		modnvec.f[k] = modn;	
	}

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
	} while( i != modnM1);

	exps.insert(exps.end(), exps.begin()+1, exps.end());
	
	for(coeffType i = 0; i < modn; i++) {
		invs[i] = exps[modn - 1 - logs[ i ] ];
	}
}
