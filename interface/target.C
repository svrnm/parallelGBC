//
// This file is part of the source of ApCoCoALib, the ApCoCoA Library.
//
//   Copyright (c) ApCoCoA Project (Prof. Dr. Martin Kreuzer, Uni Passau)
//
//   Author: 2010 Severin Neumann
//
// Visit http://apcocoa.org/ for more information regarding ApCoCoA
// and ApCoCoALib.
// Visit http://www.apcocoa.org/wiki/ApCoCoA:KnownIssues for bugs, problems 
// and known issues.
//
// ApCoCoALib is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License (version 3 or later)
// as published by the Free Software Foundation. A copy of the full
// licence may be found in the file COPYING in this directory.
//
// ApCoCoALib is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ApCoCoALib; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

#include "ApCoCoA/Algebraic/ParallelF4.H"

namespace ApCoCoA
{
namespace AlgebraicAlgorithms
{
namespace ParallelF4
{

/* SOURCE BEGIN */
<<<SOURCE>>>
/* SOURCE END */

}

unsigned long retul(RingElem r) 
{
      ZZ z;
      IsInteger(z, r);
      unsigned long n;
      IsConvertible(n, z);
      return n;
}

using namespace ParallelF4;

void ComputeParallelGroebnerBasis(const std::vector<RingElem>& gs, std::vector<RingElem>& result, int threads) {
	SparsePolyRing usedRing(AsSparsePolyRing(owner(gs[0])));
	
	// This won't work
	if (!IsQuotientRing(CoeffRing(usedRing))) {
		return;
	}

	// Identify modulus
	ideal i = DefiningIdeal(AsQuotientRing(CoeffRing(usedRing)));
	std::vector<RingElem> ge = gens(i);
	size_t modn = retul(ge[0]);
	if(modn > 32768) {
		return;
	}
	CoeffField* cf = new CoeffField( (coeffType)modn ); 

	// Get indets
	size_t n = NumIndets( AsPolyRing(usedRing) );
	TMonoid m(n);

	// Todo: Make this a choice
	TOrdering* o = new DegRevLexOrdering( n );

	vector<Polynomial> list;
	list.reserve(gs.size());
	for(size_t k = 0; k < gs.size(); k++) 
	{
		vector<Monomial> ms;
		for(SparsePolyIter spi = BeginIter(gs[k]); !IsEnded(spi); spi++)	{
			PPMonoidElem t = PP(spi);
			vector<long> values;
			exponents(values, t);
			//Coeff: retul(coeff(spi)) % modn	
			ms.push_back( make_pair( retul(coeff(spi)) % modn , m.createElement(values) ) );
		}
		Polynomial p(ms);
		//cout << p << ", ";
		list.push_back( p );
	}

        for_each(list.begin(), list.end(), bind(mem_fn(&Polynomial::order), _1, o));
        for_each(list.begin(), list.end(), bind(mem_fn(&Polynomial::bringIn), _1, cf, false));
        for_each(list.begin(), list.end(), bind(mem_fn(&Polynomial::normalize), _1, cf));

        F4 f4;

	// TODO: Make threads a choice
        vector<Polynomial> gb = f4(list, o, cf, threads);

	for(size_t k = 0; k < gb.size(); k++) {
		RingElem r(usedRing);
		for(size_t j = 0; j < gb[k].size(); j++) {
			RingElem e = RingElem(CoeffRing(usedRing), gb[k][j].first);
			PPMonoidElem m(PPM(usedRing), gb[k][j].second->getValues());
			usedRing->myPushBack( raw(r), raw(e), raw(m) );
		}
		result.push_back( r );
	}
	delete o;	
}


}
}
