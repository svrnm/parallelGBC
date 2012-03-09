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
