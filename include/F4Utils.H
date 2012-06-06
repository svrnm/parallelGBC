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
#ifndef F4_UTILS_H
#define F4_UTILS_H
#include <sys/time.h>
#include <iostream>
#include <vector>


inline double seconds()
{    
  struct timeval tv;
  gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + ((double)tv.tv_usec)/1000000.0;
}

template<typename T> std::ostream& operator<< (std::ostream& out, const std::vector<T>& v);

template<typename T> std::ostream& operator<< (std::ostream& out, const std::vector<std::vector<T> >& v);
#endif
