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
#ifndef TERM_H
#define TERM_H
#include <string.h>
#include "TMonoid.H"
#include "TOrdering.H"
#include <unordered_set>
#include <iostream>

using namespace std;

class Polynomial;

class TermInstance {
	public: 
		TMonoid* const owner;
	protected:
		degreeType degree;
		size_t hash;
		degreeType* indets;

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

		size_t getHash() const {
			return hash;
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

		const TMonoid* monoid() const {
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

		vector<long> getValues() const {
			vector<long> v;
			for(size_t i = 0; i < owner->N; i++) {
				v.push_back( indets[i] );
			}
			return v;
		}

		bool equal(const TermInstance* other) const {
			for(size_t i = 0; i < owner->N; i++) {
				if(other->indets[i] != indets[i]) return false;
			}
			return true;
		}


		friend class TMonoid;
};

class Term {
	protected:
		TMonoid* owner;
		const TermInstance* instance;
	public:
		Term(TMonoid* owner, const std::vector<degreeType>& i) : owner(owner), instance(owner->createElement(i)) { }
		Term(TMonoid* owner, const TermInstance* instance) : owner(owner), instance(instance) { }
		Term(TMonoid* owner, const std::string& s, degreeType min) : owner(owner), instance(owner->createElement(s, min)) { }

		degreeType operator[](size_t i) const {
			return instance->at(i);
		}

		degreeType deg() const {
			return instance->deg();
		}

		bool isDivisibleBy(const Term& other) const {
			return instance->isDivisibleBy(other.instance);	
		}

		Term mul(const Term& other) const {
			return Term(owner, instance->mul(other.instance));
		}

		Term lcm(const Term& other) const {
			return Term(owner, instance->lcm(other.instance));
		}

		Term div(const Term& other) const {
			return Term(owner, instance->div(other.instance));
		}

		bool equal(const Term& other) const {
			return instance->equal(other.instance);
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

		// Multiplies a set of terms with the current term
		// Uses OpenMP-Parallelization
		vector<Term> mulAll(Polynomial& in, int threads, double& timer) const;
		
		class comparator {
			public:
				const TOrdering* O;
				bool gt; // if greater

				comparator(const TOrdering* O, bool gt = false) : O(O), gt(gt) {}

				bool operator() (const Term& lhs, const Term& rhs) const
				{
					return gt ? O->cmp(rhs, lhs) < 0 : O->cmp(lhs, rhs) < 0;
				}
		};

		size_t hash() const {
			return instance->getHash();
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


ostream& operator<< (ostream &out, const Term& term);
#endif
