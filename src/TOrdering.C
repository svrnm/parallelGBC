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
#include "../include/TOrdering.H"
#include "../include/Term.H"

  int DegRevLexOrdering::cmp(const Term* a, const Term* b) const
  {   
    if(a == b) return 0;
    if(a->deg() == b->deg())
    {   
      for(size_t i = N; i > 0; i--)
      {   
        degreeType r = a->at(i-1) - b->at(i-1);
        if(r != 0)
        {   
          return r > 0 ? -1 : 1 ; 
        }   
      }   
      return 0;
    }   
    return a->deg() < b->deg() ? -1 : 1;
  }

  int LexOrdering::cmp(const Term* a, const Term* b) const 
  {
	if(a == b) return 0;
	for(long i = 0; i < N; i++)
	{
		degreeType r = a->at(i) - b->at(i);
		if(r != 0) {
			return r < 0 ? -1 : 1;
		}
	}
	return 0;
  }
