#ifndef M_ORDERING_H
#define M_ORDERING_H
#include <string>
#include "TMonoid.H"
#include "TOrdering.H"
#include "PolynomialVector.H"

using namespace std;

class MOrdering {
  public:
    const TOrdering* O;
    MOrdering(const TOrdering* O) : O(O) {}
    virtual int cmp(const TermPosition& a, const TermPosition& b) const = 0;
    virtual string name() const = 0;
};

class POTOrdering : public MOrdering
{
  public:
	POTOrdering(const TOrdering* O) : MOrdering(O) { }
  	virtual int cmp(const TermPosition& a, const TermPosition& b) const {
		if(a.second == b.second) {
			return O->cmp(a.first, b.first);
		}
		return a.second - b.second;
	}
	virtual string name() const { return "POTOrdering"; }
};
#endif
