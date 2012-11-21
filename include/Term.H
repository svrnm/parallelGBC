/**
 * This file provides the classes 'TermInstance' and 'Term'. The class 'TermInstance' is
 * only accessible via the class TMonoid (see include/TMonoid.H) or via the wrapper 'Term'.
 * A term is represented as a vector of N elements, wehre N is the number of indeterminants.
 *
 ****
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
#ifndef TERM_H
#define TERM_H
#include <stdlib.h>
#include <string>
#include "TMonoid.H"
#include "TOrdering.H"
#include <iostream>

/**
 * Forward declaration of the class Polynomial, see include/Polynomial.H
 */
class Polynomial;

class TermInstance {
	protected:
		degreeType* indets;
		TMonoid* const owner;
		size_t hash;
		degreeType degree;

		degreeType* alloc() {
			return (degreeType*)malloc(owner->N*sizeof(degreeType));
		}

		void setDegree() {
			degree = indets[0];
			for(size_t k = 1; k < owner->N; k++) {
				degree += indets[k];
			}
		}		
		void setHash() {
			hash = 0;
			for(size_t i = 0; i < owner->N; i++) {
				hash ^= indets[i] + 0x9e3779b9 + (hash << 6) + (hash >> 2);
			}
		}

		void set(size_t i, degreeType value) const {
			this->indets[i] = value;
		}

	public:
		TermInstance(TMonoid* const m) : owner(m) {
			indets = this->alloc();
		}

		TermInstance(TMonoid* const m, const degreeType* i) : owner(m)  {
			indets = alloc();
			memcpy(indets, i, owner->N*sizeof(degreeType));
			setHash();
			setDegree();
		}

		TermInstance(TMonoid* const m, const std::vector<degreeType>& v) : owner(m) {
			indets = alloc();
			for(size_t i = 0; i < owner->N; i++) {
				indets[i] = v[i];
			}
			setHash();
			setDegree();
		}

		size_t getHash() const {
			return hash;
		}

		std::vector<degreeType> getValues() const {
			std::vector<degreeType> result;
			for(size_t i = 0; i < owner->N; i++) {
				result.push_back(indets[i]);
			}
			return result;
		}

		// Multiplication constructor
		TermInstance(TMonoid* const m, const degreeType* a, const degreeType* b) : owner(m) {
			indets = alloc();
			for(size_t i = 0; i < owner->N; i++) {
				indets[i] = a[i] + b[i];
			}
			setHash();
			setDegree();
		}

		~TermInstance() {
			free(indets);
		}



		degreeType operator[](size_t i) const {
			return indets[i];
		}

		degreeType at(size_t i) const {
			return indets[i];
		}

		size_t size() const {
			return owner->N;
		}

		TMonoid* monoid() const {
			return owner;
		}

		degreeType deg() const {
			return degree;
		}

		bool isDivisibleBy(const TermInstance* other) const;

		const TermInstance* mul(const TermInstance* other) const {
			if(other->degree == 0) { return this; }
			return owner->createElement(new TermInstance(owner, indets, other->indets));
		}

		const TermInstance* mulX(size_t i) const {
			TermInstance* result = new TermInstance(owner);
			for(size_t j = 0; j < owner->N; j++) {
				result->indets[j] = indets[j];
			}
			result->indets[i]++;
			result->setDegree();
			result->setHash();
			return owner->createElement(result);
		}

		const TermInstance* lcm(const TermInstance* other) const;

		const TermInstance* div(const TermInstance* other) const {
			TermInstance* result = new TermInstance(owner);
			for(size_t i = 0; i < owner->N; i++) {
				result->indets[i] = indets[i] - other->indets[i];
			}		
			result->setDegree();
			result->setHash();	
			return owner->createElement(result);
		}

		bool equal(const TermInstance* other) const {
			for(size_t i = 0; i < owner->N; i++) {
				if(other->indets[i] != indets[i]) return false;
			}
			return true;
		}

		std::string str() const;

		friend class TMonoid;
};

class Term {
	protected:
		const TermInstance* instance;
		Term(TMonoid* owner, const TermInstance* instance) : instance(instance) { }
	public:
		Term(TMonoid* owner, const std::vector<degreeType>& i) : instance(owner->createElement(i)) { }
		Term(TMonoid* owner, const std::string& s, degreeType min) : instance(owner->createElement(s, min)) { }
		Term() {}

		degreeType operator[](size_t i) const {
			return instance->at(i);
		}

		TMonoid* monoid() const {
			return instance->monoid();
		}

		degreeType deg() const {
			return instance->deg();
		}

		bool isDivisibleBy(const Term& other) const {
			return instance->isDivisibleBy(other.instance);	
		}

		Term mul(const Term& other) const {
			return Term(instance->monoid(), instance->mul(other.instance));
		}

		Term mulX(size_t i) const {
			return Term(instance->monoid(), instance->mulX(i));
		}

		Term lcm(const Term& other) const {
			return Term(instance->monoid(), instance->lcm(other.instance));
		}

		Term div(const Term& other) const {
			return Term(instance->monoid(), instance->div(other.instance));
		}

		bool equal(const Term& other) const {
			if( instance == other.instance ) {
				return true;
			} else {
				return instance->equal(other.instance);	
			}
		}

		bool operator==(const Term& other) const {
			return instance == other.instance;
		}

		bool operator!=(const Term& other) const {
			return instance != other.instance;
		}

		size_t size() const {
			return instance->size();
		}

		class comparator {
			public:
				const TOrdering* O;
				bool gt; // if greater

				comparator() : O(NULL), gt(false) {}
				comparator(TMonoid* M, bool gt = false) : O(M->lex()), gt(gt) {}
				comparator(const TOrdering* O, bool gt = false) : O(O), gt(gt) {}

				bool operator() (const Term& lhs, const Term& rhs) const
				{
					return gt ? O->cmp(rhs, lhs) < 0 : O->cmp(lhs, rhs) < 0;
				}
		};

		size_t hash() const {
			return instance->getHash();
		}

		std::string str() const {
			return instance->str();
		}

		std::vector<degreeType> getValues() {
			return instance->getValues();
		}
};


namespace std
{    
   template <>
   struct hash<Term> : public unary_function<Term, size_t>
   {
       size_t operator()(const Term& t) const
       {
					return t.hash();
       }
   };
}


std::ostream& operator<< (std::ostream &out, const Term& term);

#endif
