#include "../include/CoeffField.H"
#include "../include/TMonoid.H"
#include "../include/Term.H"
#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <sstream>
#include <algorithm>
#include <iostream>

using namespace boost;

TMonoid::TMonoid(size_t N) : N(N), D(64/N) { 
	one = new Term(this);
	for(size_t i = 0; i < N; i++) {
		one->set(i, 0);
	}
	one->setHash();
	one->setDegree();
	terms.insert( one );
	//	terms.rehash(1024*1024*8);
}



bool TMonoid::TermEquals::operator()(const Term* const t1, const Term* const t2) const 
{
	bool r = t1->equal(t2);
	return r;
}

size_t TMonoid::TermHash::operator()(const Term* const t) const {
	return t->hash;
} 

TMonoid::~TMonoid() {
	for(TermSet::iterator it = terms.begin(); it != terms.end(); it++) {
		delete *it;
	}
	//cout << "I CREATE:\t" << iTimer << "\n";
	//cout << "F CREATE:\t" << fTimer << "\n";
}

const Term* TMonoid::createElement(Term* t)
{
	//double helper = seconds();
	pair<TermSet::iterator, bool> result = terms.insert(t);
	if(!result.second) { 
		delete t; 
		//	iTimer += seconds() - helper;
	} else {
		//	fTimer += seconds() - helper;
	}
	return *(result.first);
}

const Term* TMonoid::createElement(vector<long>& v) 
{
	degreeType* r = (degreeType*)calloc(N, sizeof(degreeType));
	long m = max(N, v.size());
	for(long i = 0; i < m; i++) {
		r[i] = (degreeType)v[i];
	}
	return createElement(new Term(this, r));
}

const Term* TMonoid::createElement(const string& s, degreeType min) { 
	degreeType* v = (degreeType*)calloc(N, sizeof(degreeType));
	if(s == "1") {
		//return createElement(new Term(this, v));
		free(v);
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
	return createElement(new Term(this, v));
}
