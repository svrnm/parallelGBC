#ifndef T_ORDERING_H
#define T_ORDERING_H
#include <string>
#include "TMonoid.H"

using namespace std;

class Term;

class TOrdering {
  public:
    const degreeType N;
    TOrdering(degreeType n) : N(n) {}
    virtual int cmp(const Term* a, const Term* b) const = 0;
    virtual string name() const = 0;
};

class DegRevLexOrdering : public TOrdering
{
  public:
	  DegRevLexOrdering(size_t N) : TOrdering(N) { }
  	virtual int cmp(const Term* a, const Term* b) const;
	  virtual string name() const { return "DegRevLex"; }
};

class LexOrdering : public TOrdering
{
  public:
	LexOrdering(size_t N) : TOrdering(N) {}
	virtual int cmp(const Term*, const Term* b) const;
	virtual string name() const { return "Lex"; }
};
#endif
