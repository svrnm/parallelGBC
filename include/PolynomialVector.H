#ifndef POLYNOMIAL_VECTOR_H
#define POLYNOMIAL_VECTOR_H
#include <iostream>
#include <vector>
#include "../include/Term.H"
#include "../include/Polynomial.H"
#include "../include/CoeffField.H"
#include "../include/TOrdering.H"

using namespace std;

typedef std::pair<const Term*, long> TermPosition;

class PolynomialVector {
	protected:
		vector<Polynomial> polynomials;

	public:
		PolynomialVector() {}
		PolynomialVector(vector<Polynomial>& ps) : polynomials(ps) { }

		Polynomial operator[](size_t i) const {
			return polynomials[i];
		}

		long size() const {
			return (long)polynomials.size();
		}

		PolynomialVector mul(const Term* t) {
			PolynomialVector result;
			for(long i = 0; i < size(); i++) {
				result.polynomials.push_back( polynomials[i].mul( t ) );
			}
			return result;
		}

		void mulBy(coeffType l, const CoeffField* field) { 
			for(long i = 0; i < size(); i++) {
				polynomials[i].mulBy(l, field);
			}	
		}

		void sub(const PolynomialVector& other, const TOrdering* O, const CoeffField* f) {
			for(long i = 0; i < size(); i++) {
				polynomials[i].sub(other.polynomials[i], O, f);
			}
		}

		// Currently only POT
		TermPosition LT() const {
			for(long i = 0; i < size(); i++) {
				if(!polynomials[i].isZero()) {
					return make_pair(polynomials[i].LT(), i);
				}
			}
			// BOOOOOOM! Ths wont work!
			return make_pair(polynomials[0].LT(),0);
		}
	
		coeffType LC() const {
			for(long i = 0; i < size(); i++) {
				if(!polynomials[i].isZero()) {
					return polynomials[i].LC();
				} 
			}
			return 0;
		}
};
#endif
