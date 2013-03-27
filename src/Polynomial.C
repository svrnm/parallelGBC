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
#include "../include/Polynomial.H"
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <boost/regex.hpp>

using namespace boost;
using namespace std;

ostream& operator<< (ostream &out, const Polynomial &poly) 
{
	bool first = true;
	// The current polynomial is zero, print 0 and return.
	if(poly.size() == 0) {
		out << 0;
		return out;
	}
	// Iterate over the monomials which are contained in poly.
	// The monomials are printed in that order which has been
	// applied to the polynomial. 
	for(size_t i = 0; i < poly.size(); i++) 
	{
		std::pair<coeffType, Term> t = poly[i];
		// Add a "+" to each polynomial, except the first
		if(!first)
		{
			out << " + ";
		}
		first = false;
		// Only print the name of the term, if it is not
		// the term of degree zero (=1)
		if(t.first != 1)
		{
			out << (int)t.first << "*";
		}
		// Print the coefficient.
		out << t.second;
	}
	return out;
}

ostream& operator<< (ostream& out, const vector<Polynomial> &polys) 
{
	out << "[";
	bool first = true;
	for(size_t i = 0; i < polys.size(); i++)
	{
		if(!first) 
		{
			out << ", ";
		}
		first = false;
		out << polys[i];
	}
	out << "]";
	return out;
}

Polynomial::Polynomial(coeffRow& cs, std::vector<Term>& ts) 
{
	if(ts.size() > 0) {
		sugarDegree = ts[0].deg();
	} else {
		sugarDegree = 0;
	}
	if(ts.size() < cs.size()) {
		cs.erase(cs.begin() + ts.size(), cs.end());
	}
	if(ts.size() > cs.size()) {
		ts.erase(ts.begin() + cs.size(), ts.end());
	}
	terms = ts;
	coeffs = cs;
}

Polynomial::Polynomial(std::vector<Monomial>& ms) 
{
	if(ms.size() > 0)
	{
		sugarDegree = ms[0].second.deg();
		for(size_t i = 0; i < ms.size(); i++) {
			coeffs.push_back(ms[i].first);
			terms.push_back(ms[i].second);
		}
	}
	else
	{
		sugarDegree = 0;
	}
}

Polynomial::Polynomial(const Term& t) 
{
	terms.push_back(t);
	coeffs.push_back(1);
	sugarDegree = t.deg();
}

Polynomial::Polynomial(std::vector<Monomial>& ms, bool purify) 
{
	if(ms.size() > 0)
	{
		sugarDegree = ms[0].second.deg();
		for(size_t i = 0; i < ms.size(); i++) {
			coeffs.push_back(ms[i].first);
			terms.push_back(ms[i].second);
		}
	}
	if(purify)
	{
		for(size_t i = 0; i < coeffs.size(); i++)
		{
			for(size_t j = i+1; j < coeffs.size();) {
				if( terms[i] == terms[j] ) {
					coeffs[i] += coeffs[j];
					coeffs.erase(coeffs.begin() + j);
					terms.erase(terms.begin() + j);
				} else {
					j++;
				}
			}
		}
	}
}


Polynomial Polynomial::mul(const Term& t) const {
	vector<Term> ts(terms.begin(), terms.end());
	coeffRow cs(coeffs.begin(), coeffs.end());
	for(size_t i = 0; i < size(); i++) {
		ts[i] = t.mul(ts[i]);
	}
	return Polynomial(cs,ts);
}

void Polynomial::mulBy(const Term& t) {
	for(size_t i = 0; i < size(); i++) {
		terms[i] = t.mul(terms[i]);
	}
}

void Polynomial::normalize(const CoeffField* field) {
	if(coeffs.size() > 0 && coeffs[0] != 0) {
		coeffType f = field->inv(coeffs[0]);
		coeffs[0] = 1;
		for(size_t i = 1; i < size(); i++) {
			coeffs[i] = field->mul(coeffs[i], f);
		}
	}
}

void Polynomial::mulBy(coeffType l, const CoeffField* field) {
	for(size_t i = 0; i < size(); i++) {
		coeffs[i] = field->mul(coeffs[i], l);
	}
}

void Polynomial::bringIn(const CoeffField* field, bool normalize) {
	for(size_t i = 0; i < coeffs.size();) {
		//coeffs[i] = ((coeffs[i] % field->modn) + field->modn) % field->modn;
		coeffs[i] = field->bringIn(coeffs[i]);
		if(coeffs[i] == 0) {
			coeffs.erase(coeffs.begin() + i);
			terms.erase(terms.begin() + i);
		} else {
			i++;
		}
	}
	if(normalize) { this->normalize(field); }
}

Polynomial Polynomial::createInstance(const string& s, TMonoid& m, degreeType min) {
	string s2 = boost::replace_all_copy(s, "-", "+-");
	boost::erase_all(s2, " ");
	vector<string> strs;
	vector<string> strs2;
	vector<Monomial> monomials;
	boost::split(strs, s2, boost::is_any_of("+"));
	for(size_t i = 0; i < strs.size(); i++) {
		if(strs[i].size() == 0) { continue; }
		string current = boost::replace_first_copy(strs[i], "x", "s");
		boost::trim(current);
		boost::split(strs2, current, boost::is_any_of("s"));
		if(strs2[0] == "") { strs2[0] = "1"; } else if(strs2[0] == "-") { strs2[0] = "-1"; }
		stringstream s1(strs2[0]);
		coeffType co;
		s1 >> co;
		string t;
		if(strs2.size() < 2) {
			t = "1";
		} else {
			strs2[1] = "x" + strs2[1];
			t = strs2[1];
		}
		monomials.push_back(make_pair(co, Term(&m, t, min)));
	}
	Polynomial result(monomials, true);
	return result;
}

// TODO: Fix this by implementing RandomAccessIterator for Polynomials
void Polynomial::order(const TOrdering* O) {
	MonomialComparator m(O);
	vector<Monomial> ms;
	for(size_t i = 0; i < terms.size(); i++) {
		ms.push_back(make_pair(coeffs[i], terms[i]));
	}
	sort(ms.begin(), ms.end(), m);
	terms.clear();
	coeffs.clear();
	for(size_t i = 0; i < ms.size(); i++) {
		terms.push_back(ms[i].second);
		coeffs.push_back(ms[i].first);
	}
}

vector<Polynomial> Polynomial::createList(const string& s, TMonoid& m, degreeType min) {
	vector<string> strs;
	boost::split(strs, s, boost::is_any_of(";,"));
	vector<Polynomial> polys;
	polys.reserve(strs.size());
	for(size_t i = 0; i < strs.size(); i++) {
		polys.push_back(Polynomial::createInstance(strs[i], m, min));
	}
	return polys;
}

bool MonomialComparator::operator()(const Monomial& lhs, const Monomial& rhs) const
{ 
	int t = O->cmp(lhs.second, rhs.second); // compare terms; 
	return t == 0 ? (lhs.first - rhs.first) > 0 : t > 0; 
}
