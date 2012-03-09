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
