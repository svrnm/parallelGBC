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
#if PGBC_USE_SSE == 1
	modnvec = __COEFF_FIELD_VECSET1( modn );
#endif

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
