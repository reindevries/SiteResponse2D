#pragma once

#define _USE_MATH_DEFINES
#define NOMINMAX

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cassert>
#include <algorithm>

//
// Debugging macros
//

#define ASSERT assert

#if _DEBUG
	#define ASSERTRETURN(x) {auto y = x; ASSERT(y); return y;}
#else
	#define ASSERTRETURN(x) {return x;}
#endif

//
// Common functions
//

template <typename T>
inline T sq(T x) {return x*x;}

template <typename T>
inline T ln(T x) {return log(x);}

template <typename T>
int sgn(T x) {return (T(0) < x) - (x < T(0));}

namespace std
{
inline size_t replace(string& str, const string& from, const string& to)
{
	size_t count = 0;
	for (size_t n = 0; (n = str.find(from, n)) != string::npos; n += to.length(), ++count)
		str.replace(n, from.length(), to);
	return count;
}

inline void to_lower(std::string& str)
	{transform(str.begin(), str.end(), str.begin(), ::tolower);}
}
