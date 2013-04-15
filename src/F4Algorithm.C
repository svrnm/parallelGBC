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
#include "../include/F4DefaultReducer.H"
#include <stdio.h>
#include <unordered_map>
#include <unordered_set>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <sstream>

using namespace std;
using namespace tbb;

namespace parallelGBC {

	void F4::updatePairs(vector<Polynomial>& polys, bool initial) 
	{
		// Setup timer to measuere how long the 'update' takes.
		double timer = F4Logger::seconds();
		// The index of the current new element in the groebner basis.
		size_t t = groebnerBasis.size();
		// Since the groebner basis grows, we store the current size
		// of the groebner basis in a variable. 
		size_t is = groebnerBasis.size();
		// Iterate over all new polynomials, which have been found in the last
		// run of reduce(...)
		for(size_t i = 0; i < polys.size(); i++)
		{
			Polynomial& h = polys[i];

			// True if the current polynomial h will be inserted into the groebner basis
			bool insertIntoG = true;

			// Check if h should be inserted: h will not be inserted into g if there
			// is already a (new!) element in the groebner basis, which has a leading
			// term that divides the leading term of h
			// Note: An old element of the groebner basis means that this element was
			// already known before the last call of reduce, so h can't have a leading
			// term which is divisible by the leading term of the old element, because
			// if this would be the case, there would have been a reduction polynomial
			// reducing this leading term.
			if(!initial) {
				for(size_t i = is; insertIntoG && i < groebnerBasis.size(); i++)
				{
					if(inGroebnerBasis[i] && h.LT().isDivisibleBy(groebnerBasis[i].LT())) 
					{
						insertIntoG = false;
					}
				}
			}

			//Check the criteria only if h will be inserted into the groebner basis
			if(insertIntoG)
			{   
				// Do the first criterium: 
				// Cancel in P all pairs (i,j) which satisfy T(i,j) = T(i,j,t), T(i,t) != T(i,j) != T(j,t) [ B_t(i,j) ]
				F4PairSet P1(pairs.key_comp());
				for(set<F4Pair>::iterator it = pairs.begin(); it != pairs.end(); it++) 
				{
					if( !it->LCM.isDivisibleBy(h.LT())  || h.lcmLT(groebnerBasis[it->i]) == it->LCM  ||  h.lcmLT(groebnerBasis[it->j]) == it->LCM  ) { 
						P1.insert( *it );
					}
				}   
				swap(pairs, P1);

				// Do the second criterium:
				// Cancel in D1 each (i,t) for which a (j,t) exists s.t. T(i,t) is a proper multiple of T(j,t) [ M(i,t) ]
				vector<bool> D1(inGroebnerBasis.begin(), inGroebnerBasis.end());
				for(size_t i = 0; i < D1.size(); i++) 
				{
					for(size_t j = i+1; D1[i] && j < D1.size(); j++)
					{
						if(D1[j])				
						{
							Term a = h.lcmLT(groebnerBasis[i]);
							Term b = h.lcmLT(groebnerBasis[j]);
							if(a != b) {
								if(a.isDivisibleBy(b))
								{
									D1[i] = false;
								}
								if(b.isDivisibleBy(a)) 
								{
									D1[j] = false;
								}
							}
						}
					}
				}

				// Do the thrd criterium:
				// In each nonvoid subset { (j,t) | T(j,t) = tau } ...
				// Attention P2 is not a multiset, so each element is unique.
				set<F4Pair, F4Pair::comparator> P2(pairs.key_comp());
				for(size_t i = 0; i < D1.size(); i++)
				{
					if(D1[i])
					{
						Term LCM = groebnerBasis[i].lcmLT(h);
						// Create a new pair: LCM, i, t, marked, sugar degree
						F4Pair newpair( LCM, i, t, LCM == groebnerBasis[i].LT().mul(h.LT()), max(groebnerBasis[i].sugar() - groebnerBasis[i].LT().deg(), h.sugar() - h.LT().deg()) + LCM.deg() );
						pair<F4PairSet::iterator,bool> ret;
						ret = P2.insert( newpair );
						// If there is a marked pair for the given LCM, store this,
						// since all pairs with this LCM will be deleted.
						if(newpair.marked && ret.second)
						{
							P2.erase(ret.first);
							P2.insert(newpair);
						}
					}
				}

				// Finally delete all (i,t) with T(i)T(j) = T(i,t).
				for(set<F4Pair, F4Pair::comparator>::iterator it = P2.begin(); it != P2.end(); it++)
				{   
					if(!it->marked)
					{ 
						pairs.insert(*it);
					}  
				}

				// Check all old elements of the groebner basis if the current element
				// divides the leading term, so the old element is reducible and can
				// removed from the result set.
				for(size_t j = 0; j < groebnerBasis.size(); j++)
				{   
					if(inGroebnerBasis[j] && groebnerBasis[j].LT().isDivisibleBy(h.LT()))
					{   
						inGroebnerBasis[j] = false;
					}
				}
				// Insert h into the groebner basis
				groebnerBasis.push_back( h );

				inGroebnerBasis.push_back( insertIntoG );
				t++;
			}
		}
		log->updateTime += F4Logger::seconds() - timer;
	}


	void F4::select() {
		std::vector<std::pair<size_t, Term> > selectedPairs;

		vector<F4Pair> tmp(pairs.begin(), pairs.end());
		
		if(withSugar) {
			sort(tmp.begin(), tmp.end(), F4Pair::sugarComparator());
			currentDegree = tmp.begin()->sugar;
		} else {
			sort(tmp.begin(), tmp.end(), F4Pair::comparator(O));
			currentDegree = tmp.begin()->LCM.deg();
		}
		if(log->verbosity & 16) {
			*(log->out) << "Degree:\t" <<currentDegree << "\n";
		}

		size_t index = 0;
		if(withSugar) {
			for(; index < tmp.size() && tmp[index].sugar == currentDegree; index++) { reducer->addSPolynomial(tmp[index].i, tmp[index].j, tmp[index].LCM); }
		} else {
			for(; index < tmp.size() && tmp[index].LCM.deg() == currentDegree; index++) { reducer->addSPolynomial(tmp[index].i, tmp[index].j, tmp[index].LCM); }
		}
		pairs.clear();
		pairs.insert(tmp.begin() + index, tmp.end());
		tmp.clear();
	}

	
	void F4::reduce(vector<Polynomial>& polys) {
		reducer->reduce(polys, currentDegree);
	}

	vector<Polynomial> F4::compute(vector<Polynomial>& generators) 
	{
		double start = F4Logger::seconds();

		pairs = F4PairSet( F4Pair::comparator(O) );

		tbb::task_scheduler_init init(threads);

		sort(generators.begin(), generators.end(), Polynomial::comparator(O, true));		

		//normalize
		for_each(generators.begin(), generators.end(), bind(mem_fn(&Polynomial::normalize), _1, field));

		updatePairs(generators, true);

		if(this->reducer == 0) {
			this->reducer = new F4DefaultReducer(this, false, 1024);
		}
		this->reducer->init();
		
		while( !pairs.empty() ) {
			vector<Polynomial> polys;
			select();
			reduce(polys);
			if(!polys.empty()) {
				updatePairs(polys, false);
			}
		}

		if(log->verbosity & 2) {
			*(log->out) << "Reduction (s): \t" << log->reductionTime << "\n";
		}
		if(log->verbosity & 4) {
			*(log->out) << "Prepare (s): \t" << log->prepareTime << "\n";
			*(log->out) << "Simplify (s): \t" << log->simplifyTime << "\n";
		}
		if(log->verbosity & 8) {
			*(log->out) << "Update (s): \t" << log->updateTime << "\n";
		}
		vector<Polynomial> result;
		for(size_t i = 0; i < groebnerBasis.size(); i++)
		{
			if(inGroebnerBasis[i])
				result.push_back(groebnerBasis[i]);
		}
		if(log->verbosity & 1) {
			*(log->out) << "Runtime (s):\t" << F4Logger::seconds() - start << "\n";
		}

		/* Enable this to collect timings
		 *(log->out) << threads << ";";
		 *(log->out) << F4Logger::seconds() - start << ";";
		 *(log->out) << reductionTime << ";";
		 *(log->out) << updateTime << ";";
		 *(log->out) << prepareTime << "\n";
		 */
		return result;
	}
}
