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
#include <sstream>
#include "../include/Term.H"
#include "../include/Polynomial.H"
#include "../include/F4Utils.H"



std::string TermInstance::str() const {
	
	std::stringstream stream;

	if(deg() > 0)
	{
		bool first = true;
		for(size_t i = 0; i < size(); i++)
		{ 
			if(indets[i] > 0)
			{ 
				if(!first)
				{ 
					stream << "*";
				}
				first = false;
				stream << "x[" << (i+1) << "]";
				if(indets[i] > 1)
				{
					stream << "^" << indets[i];
				}
			}
		}
	} else {
		stream << 1;
	}
	return stream.str();
} 

std::ostream& operator<< (std::ostream &out, const Term& term)
{ 
	out << term.str();
	return out;
}

bool TermInstance::isDivisibleBy(const TermInstance* other) const {
	if(other == this || other->degree == 0) { return true; }
	if(other->degree > degree) { return false; }
	for(size_t i = 0; i < owner->N; i++) {
		if(other->indets[i] > indets[i]) return false;
	}
	return true;
}

const TermInstance* TermInstance::lcm(const TermInstance* other) const {
	TermInstance* t = new TermInstance(this->owner);
	for(size_t i = 0; i < owner->N; i++) {
		t->indets[i] = std::max(indets[i], other->indets[i]);
	}
	t->setDegree();
	t->setHash();
	return owner->createElement(t);
}
