#ifndef T_MONOID_H
#define T_MONOID_H
#include <stdint.h>
#include <string>
#include <vector>
#include <unordered_set>
//#include <boost/unordered_set.hpp>
//#include <boost/unordered_map.hpp>
#include "F4Utils.H"

using namespace std;
//using namespace boost;

#ifndef F4_DEGREE_BITS
#ifdef F4_COEFF_BITS
#define F4_DEGREE_BITS F4_COEFF_BITS
#else
#define F4_DEGREE_BITS 32
#endif
#endif

#if F4_DEGREE_BITS <= 16
#define __MONOID_DEGREE_TYPE int16_t
#else
#define __MONOID_DEGREE_TYPE int32_t
#endif


typedef __MONOID_DEGREE_TYPE degreeType;

class Term;

class TMonoid {
	protected:
		struct TermEquals : std::binary_function<const Term* const, const Term* const , bool>
	{
		bool operator()(const Term* const t1, const Term* const t2) const ;
	};
		struct TermHash : std::unary_function<const Term* const, std::size_t>
	{
		std::size_t operator()(const Term* const t ) const;
	};

		typedef unordered_set<Term*, TermHash, TermEquals> TermSet;
		TermSet terms;
		Term* one;
		TMonoid(const TMonoid& m) : N(0), D(0) {}

	public:
		const size_t N;
		const size_t D;
		double iTimer;
		double fTimer;
		TMonoid(size_t N);
		

		const Term* createElement(const string& s, degreeType min);
		const Term* createElement(vector<long>& v);
		const Term* createElement(Term* t);

		size_t size() const {
			return terms.size();
		}
		size_t memSize() const {
			return terms.size() * ( sizeof(TMonoid*) + sizeof(degreeType)*(N+1) );
		}
		const Term* getOne() const {
			return this->one;
		}
		~TMonoid();
};
#endif
