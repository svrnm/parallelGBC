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
#include "../include/TMonoid.H"
#include "../include/Term.H"
#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <sstream>
#include <algorithm>
#include <iostream>

using namespace boost;
using namespace std;

TMonoid::TMonoid(size_t N) : N(N) { 
	one = new TermInstance(this);
	for(size_t i = 0; i < N; i++) {
		one->set(i, 0);
	}
	one->setHash();
	one->setDegree();
	//terms.insert( one );
	createElement( one );

	lexOrdering = new LexOrdering(N);
	degLexOrdering = new DegLexOrdering(N);
}



bool TMonoid::TermInstanceEquals::operator()(const TermInstance* const t1, const TermInstance* const t2) const 
{
	return (t1 == t2) || (t1 && t2 && t1->equal(t2));
}

size_t TMonoid::TermInstanceHash::operator()(const TermInstance* const t) const {
	if(t)
		return t->hash;
	else
		return 0;
} 

TMonoid::~TMonoid() {
	for(TermInstanceSet::iterator it = terms.begin(); it != terms.end(); it++) { 
		delete *it; 
	}
	delete lexOrdering;
}

const TermInstance* TMonoid::createElement(TermInstance* t)
{
	pair<TermInstanceSet::iterator, bool> result = terms.insert(t);
	if(!result.second) { 
		delete t; 
		return *(result.first);
	} else {
		return t;
	}
}

const TermInstance* TMonoid::createElement(const vector<degreeType>& v) 
{
	std::vector<degreeType> c = v;
	c.resize(N, 0);
	return createElement(new TermInstance(this, c));
}

const TermInstance* TMonoid::createElement(const string& s, degreeType min) { 
	//degreeType* v = (degreeType*)calloc(N, sizeof(degreeType));
	vector<degreeType> v(N, 0);
	if(s == "1") {
		return one;
	} else {
		std::vector<std::string> strs;
		std::vector<std::string> in_ex;
		string s2 = boost::erase_all_copy(boost::erase_all_copy(boost::erase_all_copy(s, "*"), "^"), "[");
		boost::split(strs, s2, boost::is_any_of("x"));
		for(size_t i = 1; i < strs.size(); i++) {
			boost::split(in_ex, strs[i], boost::is_any_of("]"));
			if(in_ex[1] == "") { in_ex[1] = "1"; }
			stringstream s1(in_ex[0]);
			stringstream s2(in_ex[1]);
			size_t in;
			degreeType ex;
			s1 >> in;
			s2 >> ex;
			if(in - min < N) {
				v[in-min] = ex; 
			}
		}
	}
	return createElement(new TermInstance(this, v));
}

TOrdering* TMonoid::lex() {
	return this->lexOrdering;
}

TOrdering* TMonoid::degLex() {
	return this->degLexOrdering;
}
