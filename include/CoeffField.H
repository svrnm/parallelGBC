#ifndef T_COEFF_FIELD_H
#define T_COEFF_FIELD_H
#include <assert.h>
#include <stdint.h>
#include <vector>
#include <iostream>
#ifndef F4_USE_SSE
#define F4_USE_SSE 1
#endif

#ifndef F4_COEFF_BITS 
	#define F4_COEFF_BITS 16
#endif

#if F4_COEFF_BITS <= 16
	#define __COEFF_FIELD_INTVECSIZE 8
	#define __COEFF_FIELD_GTVEC(x,y) (__builtin_ia32_pcmpgtw128 (x,y))
	#define __COEFF_FIELD_NUMBERTYPE int16_t
#else
	#define __COEFF_FIELD_INTVECSIZE 4
	#define __COEFF_FIELD_GTVEC(x,y) (__builtin_ia32_pcmpgtd128 (x,y))
	#define __COEFF_FIELD_NUMBERTYPE int32_t
#endif

using namespace std;

typedef __COEFF_FIELD_NUMBERTYPE coeffType;
// Define SSE Vector
typedef coeffType intvec __attribute__ ((vector_size (16)));
// Helper to initialize the vector
union intvecarray { intvec v; coeffType f[__COEFF_FIELD_INTVECSIZE]; };


class CoeffField /* : public CoeffField */
{ 
	public:
		const coeffType modn;

	private:
		vector<coeffType> exps;
		vector<coeffType> logs;
		vector<coeffType> invs;
		const coeffType modnM1;
		intvecarray modnvec;
		coeffType factor;
		size_t limit;
		size_t s;

	public:
		CoeffField(coeffType modn);
		coeffType mul(coeffType a, coeffType b) const
		{ 
			if( a == 0 || b == 0 ) return 0;
			return exps[logs[a] + logs[b]];
		}

		coeffType inv(coeffType a) const
		{ 
			return invs[a];
		}
		coeffType minus(coeffType a) const
		{
			if( a == 0 ) return 0;
			return modn - a;
		}
		coeffType sub(coeffType a, coeffType b) const
		{
			if( b == 0) return a;
			return b > a ? a + modn - b: a - b;
		}
		coeffType add(coeffType a, coeffType b) const
		{
			coeffType c = a + b;
			return c < modn ? c : c - modn;
		}

		coeffType getFactor(coeffType a) {
			return logs[a];
		}

		// Attention: Not safe for b equal zero !
		coeffType mulSub(coeffType a, coeffType b, coeffType c) const {
			if(b == 0) { return a; }
			b = exps[logs[b] + c];
			return b > a ? a + modn - b : a - b;
		}

		void prepare(size_t size) {
			 s = size - (size % (__COEFF_FIELD_INTVECSIZE));
			 limit = s / (__COEFF_FIELD_INTVECSIZE);
		}

		#if F4_USE_SSE
		void mulSub(vector<coeffType>& t, vector<coeffType>& o, coeffType c) const {
			vector<coeffType> tmp(o);
			if(c != 1) {
				c = logs[c];
				for(size_t k = 0; k < o.size(); k++) {
					if(o[k] != 0) { 
						tmp[k] = exps[logs[o[k]] + c]; 
					}
				}
			}
			intvec* y = (intvec*) &(tmp[0]);
			intvec* x = (intvec*) &(t[0]);
			for(size_t k = 0; k < limit; k++) {
				x[k] += (__COEFF_FIELD_GTVEC(y[k], x[k]) & modnvec.v) - y[k];
			}
			for(size_t k = s; k < t.size(); k++) {
				t[k] = tmp[k] > t[k] ? t[k] + modn - tmp[k] : t[k] - tmp[k];
			}
			tmp.clear();
		}
		#else
		void mulSub(vector<coeffType>& t, vector<coeffType>& o, coeffType c) const {
			c = logs[c];
			for(size_t k = 0; k < o.size(); k++) {
				if(o[k] != 0) {
					coeffType b = exps[logs[o[k]] + c];
					t[k] = (b > t[k]) ? t[k] - b + modn : t[k] - b;	
				}
			}
		}
		#endif
};
#endif
