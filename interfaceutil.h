#ifndef INTERF_UTIL_H
#define INTERF_UTIL_H
// ************************************************************************* //
//  File: interfaceutil.h
// ************************************************************************* //
#include <sstream>
#include <cassert>
#include<string>

using namespace std;
template <class T>
string TtoS(T arg)
{
  ostringstream os(ostringstream::out);
  os<<arg;
  return os.str(); // send to the ostringstream
}

template <class T>
T StoT(string arg)
{
  T a;

  istringstream(arg)>>a;
  return a;
}
// this function adds preceding zeros to the time such that no of digits = precision
//string precision_time(int const time,int const precision);
//extern string precision_time(int const time,int const precision);

// ************************************************************************* //
//  END: interfaceutil.h
// ************************************************************************* //
#endif
