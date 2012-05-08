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
#include "../include/F4Utils.H"

template<typename T> ostream& operator<< (ostream& out, const vector<T>& v)
{
  out << "[ ";
  for(size_t i = 0; i < v.size(); i++)
  {
    out << v[i] << " ";
  }
  out << "]";
  return out;
}

template<typename T> ostream& operator<< (ostream& out, const vector<vector<T> >& v)
{
  for(size_t i = 0; i < v.size(); i++)
  {
    out << v[i] << "\n";
  }
  return out;
}
