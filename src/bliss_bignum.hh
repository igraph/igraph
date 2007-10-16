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
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef BIGNUM_HH
#define BIGNUM_HH

#include <stdlib.h>
#include <stdio.h>
#include "bliss_defs.hh"


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
  int print(FILE *fp) {return fprintf(fp, "%Lg", v); }
};

}

#endif



#endif
