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
#ifndef T_ORDERING_H
#define T_ORDERING_H
#include <string>
#include "TMonoid.H"


class Term;

class TOrdering {
  public:
    const degreeType N;
    TOrdering(degreeType n) : N(n) {}
    virtual int cmp(const Term& a, const Term& b) const = 0;
    virtual std::string name() const = 0;
};

class DegRevLexOrdering : public TOrdering
{
  public:
	  DegRevLexOrdering(size_t N) : TOrdering(N) { }
  	virtual int cmp(const Term& a, const Term& b) const;
	  virtual std::string name() const { return "DegRevLex"; }
};

class LexOrdering : public TOrdering
{
  public:
	LexOrdering(size_t N) : TOrdering(N) {}
	virtual int cmp(const Term&, const Term& b) const;
	virtual std::string name() const { return "Lex"; }
};
#endif
