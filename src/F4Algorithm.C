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
#include "../include/F4Algorithm.H"
#include <stdio.h>
#define BREAKPOINT {while(getchar() != '\n');}

void F4::updatePairs(F4PairSet& pairs, vector<Polynomial>& polys, bool initial) 
{
	double timer = seconds();
	size_t t = groebnerBasis.size();
	for(size_t i = 0; i < polys.size(); i++)
	{
		bool insertIntoG = true;
		Polynomial& h = polys[i];
		// Cancel in P all pairs (i,j) which satisfy T(i,j) = T(i,j,t), T(i,t) != T(i,j) != T(j,t) [ B_t(i,j) ]
		F4PairSet P1(pairs.key_comp());
		for(set<F4Pair>::iterator it = pairs.begin(); it != pairs.end(); it++) 
		{
			if( !it->LCM->isDivisibleBy(h.LT())  || h.lcmLT(groebnerBasis[it->i]) == it->LCM  ||  h.lcmLT(groebnerBasis[it->j]) == it->LCM  ) { 
				P1.insert( *it );
			}   
		}   
		swap(pairs, P1);

		// Let D1 := {(i,t) | 1 <= i < t }.
		for(size_t i = 0; insertIntoG && i < groebnerBasis.size(); i++)
		{
			if(inGroebnerBasis[i] && h.LT()->isDivisibleBy(groebnerBasis[i].LT())) 
			{
				insertIntoG = false;
			}
		}

		if(insertIntoG)
		{   
			vector<bool> D1(inGroebnerBasis.begin(), inGroebnerBasis.end());
			// Cancel in D1 each (i,t) for which a (j,t) exists s.t. T(i,t) is a proper multiple of T(j,t) [ M(i,t) ]
			for(size_t i = 0; i < D1.size(); i++) 
			{
				const Term* a = h.lcmLT(groebnerBasis[i]);
				for(size_t j = 0; D1[i] && j < D1.size(); j++)
				{
					if(i != j && D1[j])				
					{
						const Term* b = h.lcmLT(groebnerBasis[j]);
						if(a->isDivisibleBy(b) && a != b)
						{
							D1[i] = false;
						}
					}
				}
			}

			// In each nonvoid subset { (j,t) | T(j,t) = tau } ...
			F4PairSet P2(pairs.key_comp());
			for(size_t i = 0; i < D1.size(); i++)
			{
				if(D1[i])
				{
					const Term* LCM = groebnerBasis[i].lcmLT(h);
					F4Pair newpair( LCM, i, t, LCM == groebnerBasis[i].LT()->mul(h.LT()), max(groebnerBasis[i].sugar() - groebnerBasis[i].LT()->deg(), h.sugar() - h.LT()->deg()) + LCM->deg() );
					pair<set<F4Pair>::iterator,bool> ret;
					// TODO: Efficient ...
					ret = P2.insert( newpair );
					if(newpair.marked && ret.second)
					{
						P2.erase(ret.first);
						P2.insert(newpair);
					}
				}
			}

			// Finally delete all (i,t) with T(i)T(j) = T(i,t).
			for(set<F4Pair>::iterator it = P2.begin(); it != P2.end(); it++)
			{   
				if(!it->marked)
				{ 
					pairs.insert(*it);
				}  
			}

			for(size_t j = 0; j < groebnerBasis.size(); j++)
			{   
				if(inGroebnerBasis[j] && groebnerBasis[j].LT()->isDivisibleBy(h.LT()))
				{   
					inGroebnerBasis[j] = false;
				}
			}
		groebnerBasis.push_back( h );
		inGroebnerBasis.push_back( insertIntoG );
		t++;
		}

	}
	//cout << "UPDATE:\t" << seconds() - timer << "\n";
	updateTime += seconds() - timer;
}

struct s
{
	bool operator() (const F4Pair& a, const F4Pair& b)
	{
		return a.sugar < b.sugar;
	}
};


void printPolyMatrix(vector<Polynomial>& v, const TOrdering* O)
{
	TermComparator t(O, true);
	set<const Term*, TermComparator> terms(t);
	for(size_t i = 0; i < v.size(); i++) {
		for(size_t j = 0; j < v[i].size(); j++) {
			terms.insert(v[i][j].second);
		}
	}


	for(set<const Term*, TermComparator>::iterator it = terms.begin(); it != terms.end(); it++)
	{
		//cout << **it << " (" << *it << ") ; ";
	}
	//cout << "\n";

	for(size_t i = 0; i < v.size(); i++) {
		set<const Term*, TermComparator>::iterator it = terms.begin();
		for(size_t j = 0; j < v[i].size(); j++) {
			for(;*it != v[i][j].second;it++) { cout << " x "; }
			cout << " " << v[i][j].first << " ";
			it++;
		}
		for(;it != terms.end(); it++) { cout << " y "; }
		cout << "\n";
	}
}

void F4::gauss(vector<vector<coeffType> >& matrix, size_t upper, vector<bool>& empty)
{
	for(size_t i = 1; i < upper; i+=2)
	{
		size_t p = 0;
		bool found = false;
		coeffType factor = 0;
		for(p = 0; !found && p < matrix[i].size(); p++) {
			found = matrix[i][p] != 0;
			factor = matrix[i][p];
		}
		p--;
		empty[i] = !found;
		if(found) {
			// Normalize
			if(factor != 1) {
				factor = field->inv(factor);
				for(size_t j = p; j < matrix[i].size(); j++) {
					matrix[i][j] = field->mul(matrix[i][j], factor);
				}
			}
			// Execute
			#pragma omp parallel for num_threads( threads )
			for(size_t j = 2; j < upper; j+=2)
			{
				size_t k = (i+j)%upper;
				if(matrix[k][p] != 0) {
					coeffType factor = field->getFactor(matrix[k][p]);
					for(size_t m = p; m < matrix[k].size(); m++)
					{
						// This is mulSub for primitives not vectors !
						matrix[k][m] = field->mulSub(matrix[k][m], matrix[i][m], factor);
					}
				}
			}
		}
	}

}

void F4::pReduce(vector<vector<F4Operation> >& ops, vector<vector<coeffType> >& rs)
{
	for(size_t i = 0; i < ops.size(); i++) {
		#pragma omp parallel for num_threads( threads ) 
		for(size_t j = 0; j < ops[i].size(); j++)
		{
			field->mulSub(rs[ops[i][j].target], rs[ops[i][j].oper], ops[i][j].factor);
		}
	}
}

size_t F4::prepare(F4PairSet& pairs, vector<Polynomial>& polys, vector<vector<F4Operation> >& ops, set<const Term*, TermComparator>& terms, vector<vector<coeffType> >& rs)
{
	double timer = seconds();
	// SELECTION
	TermComparator tog(O, true);
	vector<F4Pair> tmp(pairs.begin(), pairs.end());
	sort(tmp.begin(), tmp.end(), s());
	currentDegree = tmp.begin()->sugar;
	//std::cout << currentDegree << " degree\n";
	vector<pair<size_t, const Term*> > rows;
	size_t index;
	// Create pivots
	unordered_map<const Term*, size_t> pivots;

	for(index = 0; index < tmp.size() && tmp[index].sugar == currentDegree; index++)
	{
		rows.push_back(make_pair(tmp[index].i, tmp[index].LCM));
		rows.push_back(make_pair(tmp[index].j, tmp[index].LCM));
		pivots.insert(make_pair(tmp[index].LCM, 2*index));
	}
	pairs.clear();
	pairs.insert(tmp.begin() + index, tmp.end());
	tmp.clear();
	// SELECTION END

	//cout << index << " pairs\n";
	size_t upper = 2*index;

	unordered_map<const Term*, vector<pair<size_t, coeffType> > > pivotOps;
	unordered_set<const Term*> termsUnordered;

	vector<vector<Monomial> > rightSide;
	rightSide.reserve(rows.size());

	double testtimer = 0;

	for(size_t i = 0; i < rows.size(); i++) 
	{
		size_t currentRow = rows[i].first;
		rightSide.push_back( vector<Monomial>() );
		//rows[i].second->divToVector(groebnerBasis[currentRow].LT(), ir);
		// For the pivot rows (even rows and lower part) we start at 1 

		// precalculated monomials
		// 50%
		const Term* ir = rows[i].second->div(groebnerBasis[currentRow].LT());
		vector<const Term*> pcm = ir->mulAll(groebnerBasis[currentRow], threads, testtimer);

		// 30%	
		for(size_t j =  (i > upper || i % 2 == 0 ? 1 : 0);  j < groebnerBasis[currentRow].size() ; j++) 
		{
			coeffType coeff = groebnerBasis[currentRow][j].first;
			const Term* t = pcm[j];
			bool wontFound = termsUnordered.count(t) > 0;	
			bool found = false;

			// If there is not yet a pivot for t, try to create one
			if(!wontFound) {
				found = pivots.count(t) > 0;
				if(!found)
				{
					for(size_t k = 0; !found && k < groebnerBasis.size(); k++) 
					{
						if(inGroebnerBasis[k] && t->isDivisibleBy(groebnerBasis[k].LT())) {
							found = true;
							rows.push_back(make_pair(k, t));
							pivots.insert(make_pair(t, rows.size()-1));
						}
					}
				}
			}

			// Eliminate if possible
			if(found) {
				pivotOps[t].push_back( make_pair(i, coeff) );
			} else {
				rightSide[i].push_back( make_pair(coeff, t) );
				if(!wontFound) termsUnordered.insert(t);
			}
		}
	}

	terms.insert(termsUnordered.begin(), termsUnordered.end());

	//rs.assign(rightSide.size(), vector<coeffType>(terms.size(), 0) );

	size_t pad = __COEFF_FIELD_INTVECSIZE;
	rs.assign(rightSide.size(), vector<coeffType>( (( terms.size()+pad-1 )/ pad ) * pad, 0) );
	for(size_t i = 0; i < rightSide.size(); i++) {
		size_t j = 0;
		size_t k = 0;
		for(set<const Term*, TermComparator>::iterator it = terms.begin(); j < rightSide[i].size() /*&& it != terms.end()*/; it++) {
			if(rightSide[i][j].second == *it) {
				rs[i][k] = rightSide[i][j].first ;
				j++;
			}
			k++;
		}
	}
	rightSide.clear();
	
	map<const Term*, vector<pair<size_t, coeffType> >, TermComparator> pivotOpsOrdered(pivotOps.begin(), pivotOps.end(), tog);
	pivotOps.clear();
	ops.push_back( vector<F4Operation >() );

	#if 1 
	vector<size_t> l(rs.size(),0);
	for(map<const Term*, vector<pair<size_t, coeffType> >, TermComparator>::reverse_iterator it = pivotOpsOrdered.rbegin(); it != pivotOpsOrdered.rend(); it++)
	{
		size_t o = pivots[it->first];
		for(size_t i = 0; i < it->second.size(); i++)
		{
			size_t t = it->second[i].first;
			if(l[ o ] > l[t]) {
				l[t] = l[o];
			}
			ops[ l[t] ].push_back( F4Operation(t,o,it->second[i].second) );
			l[t]++; // one operation per level, attention this also affects the following if-statements

			if(l[t] >= ops.size()) {
				ops.push_back( vector<F4Operation>() );
			}
		}
	}
	#else
	size_t l = 0;
	for(map<const Term*, vector<pair<size_t, coeffType> >, TermComparator>::reverse_iterator it = pivotOpsOrdered.rbegin(); it != pivotOpsOrdered.rend(); it++)
	{
		size_t o = pivots[it->first];
		for(size_t i = 0; i < it->second.size(); i++)
		{
			size_t t = it->second[i].first;
			ops[ l ].push_back( F4Operation(t, o, it->second[i].second) );
		}
		l++;
		ops.push_back( vector<F4Operation>() );
	}
	#endif
	pivotOpsOrdered.clear();
	ops.pop_back();
	prepareTime += seconds() - timer;
	//cout << "Operations:\t" << ops.size() << "\n";
	//cout << "Preparation:\t" << seconds() - timer << "\n";
	return upper;
}

void F4::reduce(F4PairSet& pairs, vector<Polynomial>& polys)
{
	TermComparator tog(O, true);
	set<const Term*, TermComparator> terms(tog);
	vector<vector<F4Operation> > ops;
	vector<vector<coeffType> > rs;
	//vector<vector<size_t> > setOffset; 
	size_t upper = prepare(pairs, polys, ops, terms, rs);

	// ELIMINATE
	double timer = seconds();
	pReduce(ops, rs);
	ops.clear();
	//timer = seconds();

	vector<bool> empty(upper, false); // too large, FIX?
	
	gauss(rs, upper, empty);
	//cout << "GAUSSR:\t" << seconds()-timer << "\n";
	// ELIMINATE END
	reductionTime += seconds()-timer;
	//cout << "REDUCE:\t" << seconds()-timer << "\n";
	for(size_t i = 1; i < upper; i+=2)
	{
		if(!empty[i])
		{
			Polynomial p(currentDegree);
			size_t j = 0;
			for(set<const Term*, TermComparator>::iterator it = terms.begin(); it != terms.end(); it++) 
			{
				if(rs[i][j] != 0)
				{
					p.push_back(make_pair(rs[i][j], *it));
				}
				j++;
			}
			polys.push_back( p );
		}
	}
	//cout << polys.size() << " new elements\n";
	//postReduce(polys);
}

void F4::postReduce(vector<Polynomial>& polys) 
{

}


vector<Polynomial> F4::operator()(vector<Polynomial>& generators, const TOrdering* o, CoeffField* field, int threads)
{
	double start = seconds();
	this->field = field;
	this->threads = threads;
	this->O = o;
	F4PairComparator f4pc(o);
	F4PairSet pairs(f4pc);
	updateTime = 0;
	prepareTime = 0;
	reductionTime = 0;

	sort(generators.begin(), generators.end(), PolynomialComparator(o, true));		
	//normalize
	for_each(generators.begin(), generators.end(), bind(mem_fn(&Polynomial::normalize), _1, field));

	updatePairs(pairs, generators,true);

	while( !pairs.empty() ) {
		vector<Polynomial> polys;
		reduce(pairs, polys);
		if(!polys.empty()) {
			updatePairs(pairs, polys);
		}
		//cout << pairs.size() << " pairs\n";
		/*size_t t = 0;	
		size_t k = 0;
		for(size_t i = 0; i < groebnerBasis.size(); i++)
		{
			t += groebnerBasis[i].size();	
			if(inGroebnerBasis[i]) {
				k++;
			}
		}*/
	}
	cout << "Reduction (s): \t" << reductionTime << "\n";
	//cout << "P-Time: \t" << prepareTime << "\n";
	//cout << "U-Time: \t" << updateTime << "\n";
	vector<Polynomial> result;
	for(size_t i = 0; i < groebnerBasis.size(); i++)
	{
		if(inGroebnerBasis[i])
			result.push_back(groebnerBasis[i]);
	}
	cout << "Runtime (s):\t" << seconds() - start << "\n";
	return result;
}
