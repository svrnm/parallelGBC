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
	//cout << exps[13] << "\n";
	//cout << exps[13+modnM1] << "\n";

	for(coeffType i = 0; i < modn; i++) {
		invs[i] = exps[modn - 1 - logs[ i ] ];
	}
}
