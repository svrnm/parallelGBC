#include "../include/Polynomial.H"
//#include <algorithm>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <boost/regex.hpp>

using namespace boost;

ostream& operator<< (ostream &out, const Polynomial &poly) 
{
	bool first = true;
	if(poly.size() == 0) {
		out << 0;
		return out;
	}
	for(size_t i = 0; i < poly.size(); i++) 
	{
		std::pair<coeffType, const Term*> t = poly[i];
		if(!first)
		{
			out << " + ";
		}
		first = false;
		if(t.first != 1)
		{
			out << t.first << "*";
		}
		out << *t.second;
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

Polynomial Polynomial::mul(const Term* t) const {
	vector<Monomial> ms(monomials.begin(), monomials.end());
	for(size_t i = 0; i < size(); i++) {
		//ts.push_back(t->mul(terms[i]));
		ms[i].second = t->mul(ms[i].second);
	}
	return Polynomial(ms);
}

void Polynomial::mulBy(const Term* t) {
	for(size_t i = 0; i < size(); i++) {
		monomials[i].second = t->mul(monomials[i].second);
	}
}

void Polynomial::normalize(const CoeffField* field) {
	if(monomials.size() > 0 && monomials[0].first != 0) {
		coeffType f = field->inv(monomials[0].first);
		monomials[0].first = 1;
		for(size_t i = 1; i < size(); i++) {
			monomials[i].first = field->mul(monomials[i].first, f);
		}
	}
}

void Polynomial::mulBy(coeffType l, const CoeffField* field) {
	for(size_t i = 0; i < size(); i++) {
		monomials[i].first = field->mul(monomials[i].first, l);
	}
}

void Polynomial::bringIn(const CoeffField* field, bool normalize) {
	for(size_t i = 0; i < monomials.size();) {
		monomials[i].first = ((monomials[i].first % field->modn) + field->modn) % field->modn;
		if(monomials[i].first == 0) {
			monomials.erase(monomials.begin() + i);
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
		monomials.push_back(make_pair(co, m.createElement(t, min)));
	}
	Polynomial result(monomials, true);
	return result;
}

void Polynomial::order(const TOrdering* O) {
	MonomialComparator m(O);
	sort(monomials.begin(), monomials.end(), m);
}

void Polynomial::sub(const Polynomial& other, const TOrdering* O, const CoeffField* f) {
	vector<Monomial>::iterator m1 = monomials.begin();
	vector<Monomial>::const_iterator m2 = other.monomials.begin();
	vector<Monomial> result;
	back_insert_iterator<vector<Monomial> > r(result);

	while(m1 != monomials.end() && m2 != other.monomials.end()) {
		int c = O->cmp(m1->second, m2->second);
		if(c == 0) {
			coeffType coeff = f->sub(m1->first, m2->first);
			if(coeff != 0) {
				r = make_pair(coeff, m1->second);
			}
			m1++;
			m2++;
		} else if (c > 0) {
			r = *m1;
			m1++;
		} else {
			r = make_pair(f->minus(m2->first), m2->second);
			m2++;
		}
	}
	// Either m1 is end or m2, never both!
	if(m1 == monomials.end()) {
		for(; m2 != other.monomials.end(); m2++) {
			r = make_pair(f->minus(m2->first), m2->second);
		}
	} else if(m2 == other.monomials.end()) {
		copy(m1,monomials.end(),r);
	}
	swap(monomials, result);
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
