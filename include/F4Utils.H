#ifndef F4_UTILS_H
#define F4_UTILS_H
#include <sys/time.h>
#include <iostream>
#include <vector>

using namespace std;

inline double seconds()
{    
  struct timeval tv;
  gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + ((double)tv.tv_usec)/1000000.0;
}

template<typename T> ostream& operator<< (ostream& out, const vector<T>& v);

template<typename T> ostream& operator<< (ostream& out, const vector<vector<T> >& v);
#endif
