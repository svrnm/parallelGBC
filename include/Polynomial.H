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
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H
#include <iostream>
#include <vector>
#include "../include/Term.H"
#include "../include/CoeffField.H"
#include "../include/TOrdering.H"

typedef std::pair<coeffType, Term> Monomial;

class Polynomial {
	private:
		coeffRow coeffs;
		std::vector<Term> terms;
		degreeType sugarDegree;

	public:
		Polynomial(degreeType s) : sugarDegree(s) { }

		// Create empty polynomial (equal to 0)
		Polynomial() {
			sugarDegree = 0;
		}
		Polynomial(coeffRow& cs, std::vector<Term>& ts);
		Polynomial(std::vector<Monomial>& ms);
		Polynomial(const Term& t);
		Polynomial(std::vector<Monomial>& ms, bool purify);

		static Polynomial createInstance(const std::string& s, TMonoid& m, degreeType min = 1);
		static std::vector<Polynomial> createList(const std::string& s, TMonoid& m, degreeType min = 1);

		Monomial operator[](size_t i) const {
			return std::make_pair(coeffs[i], terms[i]);
		}

		coeffType coeff(size_t i) const {
			return coeffs[i];
		}

		Term term(size_t i) const {
			return terms[i];
		}

		size_t size() const {
			return coeffs.size();
		}

		size_t capacity() const {
			return terms.capacity();
		}

		Polynomial mul(const Term& t) const;

		void mulBy(const Term& t);

		bool isZero() const {
			return coeffs.size() < 1 || coeffs[0] == 0; // Attention 0*LT+k*T with k != 0 will break this!
		}

		void normalize(const CoeffField* field);

		void mulBy(coeffType l, const CoeffField* f);

		void bringIn(const CoeffField* field, bool normalize = true);

		void order(const TOrdering* O);

		Term lcmLT(const Polynomial& other) const {
			return LT().lcm(other.LT());
		}

		Term LT() const {
			return terms[0];
		}

		coeffType LC() const {
			return coeffs[0];
		}

		void setSugar(degreeType s) {
			sugarDegree = s;
		}

		degreeType sugar() const {
			return sugarDegree;
		}

		void push_back(const Monomial& m) {
			coeffs.push_back(m.first);
			terms.push_back(m.second);
		}

		void push_back(coeffType c, const Term& t) {
			coeffs.push_back(c);
			terms.push_back(t);
		}

		std::vector<Term> support() {
			return terms;
		}

		struct comparator {
			const TOrdering* O;
			bool gt; // if greater

			comparator(const TOrdering* O, bool gt = false): O(O), gt(gt) {}

			bool operator() (const Polynomial& lhs, const Polynomial& rhs) const 		{
				return gt ? O->cmp(rhs.LT(), lhs.LT()) < 0 : O->cmp(lhs.LT(), rhs.LT()) < 0;
			}
		};
};

std::ostream& operator<< (std::ostream &out, const Polynomial &poly);

class MonomialComparator {
	public:
		const TOrdering* O;
		std::greater<coeffType> g;


		MonomialComparator(const TOrdering* O) : O(O) {}

		bool operator() (const Monomial& lhs, const Monomial& rhs) const;
};

#endif
