/*
 Copyright (C) 2003-2006 Tommi Junttila

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License version 2
 as published by the Free Software Foundation.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/* FSF address fixed in the above notice on 1 Oct 2009 by Tamas Nepusz */

#ifndef BLISS_BIGNUM_HH
#define BLISS_BIGNUM_HH

using namespace std;

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <sstream>
#include "bliss_defs.hh"
#include "igraph_math.h"

#include "memory.h"
#include "error.h"

/*
 * Simple class for big integers (or approximation of such) in order
 * compute group sizes.
 * Set BLISS_USE_GMP in defs.hh to use the GMP library.
 */


#if defined(BLISS_USE_GMP)

#include <gmp.h>

namespace igraph {

class BigNum
{
  mpz_t v;
public:
  BigNum() {mpz_init(v); }
  ~BigNum() {mpz_clear(v); }
  void assign(const int n) {mpz_set_si(v, n); }
  void multiply(const int n) {mpz_mul_si(v, v, n); }
  int print(FILE *fp) {return mpz_out_str(fp, 10, v); }
  int tostring(char **str) { 
    *str=igraph_Calloc(mpz_sizeinbase(v, 10)+2, char);
    if (! *str) { 
      IGRAPH_ERROR("Cannot convert big number to string", IGRAPH_ENOMEM);
    }
    mpz_get_str(*str, 10, v);
    return 0;
  }
};

}

#else

namespace igraph {

class BigNum
{
  long double v;
public:
  BigNum(): v(0.0) {}
  void assign(const int n) {v = (long double)n; }
  void multiply(const int n) {v *= (long double)n; }
#if defined WIN32 || defined WIN64 || defined _Win32 || defined _Win64
  int print(FILE *fp) {return fprintf(fp, "%g", (double)v); }
#else
  int print(FILE *fp) {return fprintf(fp, "%Lg", v); }
#endif
  int tostring(char **str) {
    int size=static_cast<int>( (log(abs(v))/log(10.0))+4 );
    *str=igraph_Calloc(size, char );
    if (! *str) {
      IGRAPH_ERROR("Cannot convert big number to string", IGRAPH_ENOMEM);
    }
    stringstream ss;
    ss << v;
    strncpy(*str, ss.str().c_str(), size);
    return 0;
  }
};

}

#endif



#endif
