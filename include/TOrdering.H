/**
 * This file provides the abstract type of a term ordering and implementations of this
 * type. Currently there are the degree reverse lexicographic term ordering, the lexicographic
 * term ordering and the degree lexicographic term ordering.
 *
 * A given term ordering O gives for two terms a and b the information if a is less than, equal
 * or greater than b by calling the function cmp(). The result is an integer, where 0 means
 * equal, a positive value means graeter than and a negative value less than.
 ******
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
#ifndef T_ORDERING_H
#define T_ORDERING_H
#include <string>
#include "TMonoid.H"

/**
 * Forward declaration of the class Term. See include/Term.H for details
 */
class Term;

/**
 * This class represents an abstract term ordering. A term ordering is given
 * for a certain number of indeterminants (N) and provides a compare function
 * cmp() which gives the information if a term a is equal, greater than or
 * less than a term b.
 */
class TOrdering {
  public:
		/**
		 * The number of indeterminants which should be compared.
		 * Attention, if one of the given terms has less than N
		 * indeterminants the comparision will misbehave.
		 */
    const degreeType N;
		
		/**
		 * Construct a new instance of term ordering.
		 */
    TOrdering(degreeType n) : N(n) {}

		/**
		 * Compare the two terms a and b using the given implementation.
		 *
		 * If a<b return a negative value
		 * If a>b return a positive value
		 * If a=b return 0.
		 */
    virtual int cmp(const Term& a, const Term& b) const = 0;
		
		/**
		 * Return the name of the given term ordering
		 */
    virtual std::string name() const = 0;
};

/**
 * This class represents the degree reverse lexicographic term ordering. A
 * term a is less than a term b, if 
 * - the degree of a is less than the degree of b
 * - or if the degree is equal the indeterminants are compared in reverse lexicographic ordering,
 *   which means that a is less than b if a has a larger value in an indeterminant close to the
 *   end than b.
 */
class DegRevLexOrdering : public TOrdering
{
  public:
		/**
		 * @inherit
		 */
	  DegRevLexOrdering(size_t N) : TOrdering(N) { }

		/**
		 * @inherit
		 */
  	virtual int cmp(const Term& a, const Term& b) const;

		/**
		 * @inherit
		 */
	  virtual std::string name() const { return "DegRevLex"; }
};

/**
 * This class represents the lexciographic term ordering. A term a is less than a term b, if
 * a has a smaller value in an indeterminant close to the beginning than b.
 */
class LexOrdering : public TOrdering
{
  public:
	/**
	 * @inherit
	 */
	LexOrdering(size_t N) : TOrdering(N) {}

	/**
	 * @inherit
	 */
	virtual int cmp(const Term&, const Term& b) const;

	/**
	 * @inherit
	 */
	virtual std::string name() const { return "Lex"; }
};

/**
 * This class represents the degree lexicographic term ordering. A term a is less than a term b
 * if 
 * - the degree of a is less than the degree of b
 * - or the degree is equal and a is less than b regarding the lexicographic term ordering.
 */
class DegLexOrdering : public TOrdering
{
	public:
	/**
	 * @inherit
	 */
	DegLexOrdering(size_t N) : TOrdering(N) {}

	/**
	 * @inherit
	 */
	virtual int cmp(const Term&, const Term& b) const;

	/**
	 * @inherit
	 */
	virtual std::string name() const { return "DegLex"; }
};
#endif
