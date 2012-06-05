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
#include "../include/Term.H"
#include "../include/Polynomial.H"
#include "../include/F4Utils.H"

ostream& operator<< (std::ostream &out, const Term& term)
{ 
	if(term.deg() > 0)
	{
		bool first = true;
		for(size_t i = 0; i < term.size(); i++)
		{ 
			if(term[i] > 0)
			{ 
				if(!first)
				{ 
					out << "*";
				}
				first = false;
				out << "x[" << (i+1) << "]";
				if(term[i] > 1)
				{ 
					out << "^" << term[i];
				}
			}
		}
	} else {
		out << 1;
	}
	return out;
} 

bool TermInstance::isDivisibleBy(const TermInstance* other) const {
	if(other == this) { return true; }
	if(other->degree > degree) { return false; }
	for(size_t i = 0; i < owner->N; i++) {
		if(other->indets[i] > indets[i]) return false;
	}
	return true;
}

const TermInstance* TermInstance::lcm(const TermInstance* other) const {
	TermInstance* t = new TermInstance(this->owner);
	for(size_t i = 0; i < owner->N; i++) {
		t->indets[i] = max(indets[i], other->indets[i]);
	}
	t->setDegree();
	t->setHash();
	return owner->createElement(t);
}

vector<Term> Term::mulAll(Polynomial& in, int threads, double& timer) const {
	vector<Term> result;

	if(this->deg() == 0) {
		for(size_t i = 0; i < in.size(); i++) {
			result.push_back( in[i].second );
		}
		return result;
	}	

	//double helper = seconds();
	//vector<TermInstance*> tmp(in.size(), NULL);
	// This doesn't really speedup (the critical section takes the most computation time)
	//#pragma omp parallel for num_threads( 1 ) 
	for(size_t i = 0; i < in.size(); i++) {
		//result[i] = Term(owner, in[i].second->indets, indets);
		result.push_back( in[i].second.mul( *this ) );
		//#pragma omp critical
		//result[i] = owner->createElement(tmp[i]);
	}
	//timer += seconds() - helper;
	return result;
}
