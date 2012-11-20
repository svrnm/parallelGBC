/**
 * Header fiel for the polynomials. Polynomials are implemented as holder of two vectors
 * containing terms and their matching coefficients. The polynomials do not know anything
 * about the coefficient field and the term ordering. So the functions order(...) and bringIn(...)
 * are required to put a polynomial into a given ordering or field. This allwos to switch
 * the ordering and the field at runtime, however not the TermMonoid, i.e. the number of indeterminants.*
 *
 * --
 * * If you look at the implementation it is possibile to mix the terms of different monoids, which is
 * of course a bad idea...
 *
 ***************
 *
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

/**
 * Monomials are a pair of a coefficient and a term, e.g.:
 *
 * 3*(x[1]^2*x[2]^3)
 */
typedef std::pair<coeffType, Term> Monomial;

/**
 * This class represents a polynomial. A polynomial is defined by the terms
 * and the coefficients. Furthermore for groebner basis computation the sugar
 * degree (http://dl.acm.org/citation.cfm?id=120701) of a polynomial can be set.
 * The logic of the sugar degree is done mostly within the groebner basis computation.
 *
 */
class Polynomial {
	private:
		/**
		 * The coefficients of the elements in the polynomial. The coefficient at position
		 * i is the coefficient of the term at position i.
		 */
		coeffRow coeffs;
		/**
		 * The terms of the elements in the polynomial. The term at position i has the
		 * coefficient at position i.
		 */
		std::vector<Term> terms;
		/**
		 * The sugar degree of the given polynomial.
		 */
		degreeType sugarDegree;

	public:
		/**
		 * Constructs a empty polynomial with sugar degree s. The elements of the polynomial
		 * can be pushed afterwards using the push_back(...) method.
		 */
		Polynomial(degreeType s) : sugarDegree(s) { }

		/**
		 * Construct an empty polynomial with sugar degree zero (which is equal to the 0-polynomial)
		 */
		Polynomial() {
			sugarDegree = 0;
		}
		/**
		 * Construct a polynomial which has coefficients cs and terms ts
		 */
		Polynomial(coeffRow& cs, std::vector<Term>& ts);
		/**
		 * Construct a polynomial whit monomials ms. The vector of monomials
		 * will be decomposed into terms and coefficients.
		 */
		Polynomial(std::vector<Monomial>& ms);
		/**
		 * Construct a polynomial which is 1*t, where t is a term
		 */
		Polynomial(const Term& t);
		/**
		 * Construct a polynomial with monomials ms. The vector of monomials
		 * will be decomposed into terms and coefficients. If the option 'purify'
		 * is enabled, the polynomial will be checked for not-unique terms.
		 */
		Polynomial(std::vector<Monomial>& ms, bool purify);

		/**
		 * Returns a polynomial which is constructed by parsing the given string s.
		 * The terms will be elements of the Monoid m. The parameter 'min' presets
		 * the smallest x[i] within s. 
		 *
		 * The given string 's' should have the following form:
		 * - Indeterminats are given as x[i], where is is an positive integer
		 * - Powers of indeterminats are given as x[i]^n where n is an positive integer
		 * - Multiple indeterminats are multiplied by '*': x[i]^n*x[j]^m*...
		 * - Coefficients are given as positive integers and are multiplied to a term, e.g. 5*x[1]^3*x[2]^4
		 * - Multiple coefficients and terms are joined by '+' or '-', e.g. 5*x[1]^3*x[2]^4+8*x[3]^4-4*x[1]+1
		 */
		static Polynomial createInstance(const std::string& s, TMonoid& m, degreeType min = 1);

		/**
		 * Returns a list of polynomials which is constructed by parsing the given string s.
		 * The terms will be elements of the Monoid m. The parameter 'min' presets the smallest
		 * x[i] within s.
		 *
		 * The given string 's' should have the following form:
		 * - polynomials are represented as described above (createInstance)
		 * - polynomials are separeted by ", ": p0, p1, p2, ...
		 */
		static std::vector<Polynomial> createList(const std::string& s, TMonoid& m, degreeType min = 1);

		/**
		 * Return the monomial at position i
		 */
		Monomial operator[](size_t i) const {
			return std::make_pair(coeffs[i], terms[i]);
		}

		/**
		 * Return the coefficient at position i
		 */
		coeffType coeff(size_t i) const {
			return coeffs[i];
		}

		/**
		 * Return the term at position i
		 */
		Term term(size_t i) const {
			return terms[i];
		}

		/**
		 * Return the 'size' of the polynomial which is equvialent to the number of terms+coefficients contained 
		 * in the polynomial
		 */
		size_t size() const {
			return coeffs.size();
		}

		/**
		 * Multiply all terms of the polynomial with the given term 't' and return a new polynomial.
		 * The result is for a polynomial 'p' the new polynomial 'p*t'.
		 */
		Polynomial mul(const Term& t) const;

		/**
		 * Multply all terms of the polynomial with the given term 't'. The given polynomial is
		 * changed to 'p*t'
		 */
		void mulBy(const Term& t);

		/**
		 * Return true if the given polynomial is zero.
		 */
		bool isZero() const {
			return coeffs.size() < 1 || coeffs[0] == 0; // Attention 0*LT+k*T with k != 0 will break this!
		}

		/**
		 * Normalize the polynomial using the given field. This means that the current leading coefficient
		 * (coeffs[0]) has the value 1 and all other coefficients are multiplied with the inverse
		 * of the old leading coefficient
		 */
		void normalize(const CoeffField* field);

		/**
		 * Multiply all coefficients by the given coefficient 'l' within the field 'f'. Attention this
		 * requires the values of the coefficients to be in 'f'. If this is not the case, the computation
		 * may misbehave.
		 */
		void mulBy(coeffType l, const CoeffField* f);

		/**
		 * Bring all coefficients into the 'field'. Set the 'normalize' to true, to normalize the coefficients
		 * afterwards.
		 */
		void bringIn(const CoeffField* field, bool normalize = true);

		/**
		 * Apply the term ordering 'O' to this polynomial. This will rearrange the terms (and coefficients).
		 * At the end the term at the first position (terms[0]) is leading and all other terms are stored
		 * in deacreasing order.
		 */
		void order(const TOrdering* O);

		/**
		 * Return the least common multiple of the leading term of this and the leading term of other
		 */
		Term lcmLT(const Polynomial& other) const {
			return LT().lcm(other.LT());
		}

		/**
		 * Return the leading term
		 */
		Term LT() const {
			return terms[0];
		}

		/**
		 * Return the leading coefficient
		 */
		coeffType LC() const {
			return coeffs[0];
		}

		/**
		 * Set the sugar degree of the current polynomial
		 */
		void setSugar(degreeType s) {
			sugarDegree = s;
		}

		/**
		 * Get the sugar degree of the current polynomial
		 */
		degreeType sugar() const {
			return sugarDegree;
		}

		/**
		 * Add a the monomial m to the end of the polynomial.
		 * Ordering of the terms will not be reorganized
		 */
		void push_back(const Monomial& m) {
			coeffs.push_back(m.first);
			terms.push_back(m.second);
		}

		/**
		 * Add a pair of a term and a coefficient to the end of the polynomial.
		 * Ordering of the terms will not be reorganized.
		 */
		void push_back(coeffType c, const Term& t) {
			coeffs.push_back(c);
			terms.push_back(t);
		}

		/**
		 * Return a copy of all terms of the polynomial
		 */
		std::vector<Term> support() {
			return terms;
		}

		/**
		 * A polynomial comparator. The comparator compares the leading terms of two
		 * given polynomials using the given term ordering 'O'. The flag 'gt' can be
		 * set to change the comparator from "less than" to "greater than".
		 *
		 * Of course this ordering, is not complete, since two polynomials having the
		 * same leading term may be different in their tail
		 */
		struct comparator {
			/**
			 * The ordering which is used for the comparison
			 */
			const TOrdering* O;
			/**
			 * True if this is a "greater than " comparator
			 */
			bool gt;

			/**
			 * Construct a new comparator object.
			 */
			comparator(const TOrdering* O, bool gt = false): O(O), gt(gt) {}

			/**
			 * Return true if the leading term of lhs is less than the leading term
			 * of rhs. if gt is set to true, true is returned if lhs is greater than
			 * then rhs.
			 */
			bool operator() (const Polynomial& lhs, const Polynomial& rhs) const 		{
				return gt ? O->cmp(rhs.LT(), lhs.LT()) < 0 : O->cmp(lhs.LT(), rhs.LT()) < 0;
			}
		};
};

/**
 * Print a polynomial to a output stream
 */
std::ostream& operator<< (std::ostream &out, const Polynomial &poly);

/**
 * A monomial comparator. The comparator compres the leading terms of two
 * monomials using the given term ordering 'O'.
 */
class MonomialComparator {
	public:
		/**
		 * The ordering which is used for comparison
		 */
		const TOrdering* O;

		/**
		 * Construct a new monomial comparator object.
		 */
		MonomialComparator(const TOrdering* O) : O(O) {}

		/**
		 * Return true if the leading term of lhs is less than the leading
		 * term of rhs
		 */
		bool operator() (const Monomial& lhs, const Monomial& rhs) const;
};

#endif
