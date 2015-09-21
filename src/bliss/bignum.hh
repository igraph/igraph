#ifndef BLISS_BIGNUM_HH
#define BLISS_BIGNUM_HH

/*
  Copyright (c) 2003-2015 Tommi Junttila
  Released under the GNU Lesser General Public License version 3.
  
  This file is part of bliss.
  
  bliss is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, version 3 of the License.

  bliss is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with bliss.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <sstream>
#include "defs.hh"

#include "igraph_memory.h"
#include "igraph_error.h"

#if defined(BLISS_USE_GMP)
#include <gmp.h>
#endif

namespace bliss {

/**
 * \brief A very simple class for big integers (or approximation of them).
 *
 * If the compile time flag BLISS_USE_GMP is set,
 * then the GNU Multiple Precision Arithmetic library (GMP) is used to
 * obtain arbitrary precision, otherwise "long double" is used to
 * approximate big integers.
 */

#if defined(BLISS_USE_GMP)

class BigNum
{
  mpz_t v;
public:
  /**
   * Create a new big number and set it to zero.
   */
  BigNum() {mpz_init(v); }

  /**
   * Destroy the number.
   */
  ~BigNum() {mpz_clear(v); }

  /**
   * Set the number to \a n.
   */
  void assign(const int n) {mpz_set_si(v, n); }

  /**
   * Multiply the number with \a n.
   */
  void multiply(const int n) {mpz_mul_si(v, v, n); }

  /**
   * Print the number in the file stream \a fp.
   */
  size_t print(FILE* const fp) const {return mpz_out_str(fp, 10, v); }

  int tostring(char **str) const {
    *str=igraph_Calloc(mpz_sizeinbase(v, 10)+2, char);
    if (! *str) {
      IGRAPH_ERROR("Cannot convert big number to string", IGRAPH_ENOMEM);
    }
    mpz_get_str(*str, 10, v);
    return 0;
  }

};

#else

class BigNum
{
  long double v;
public:
  /**
   * Create a new big number and set it to zero.
   */
  BigNum(): v(0.0) {}

  /**
   * Set the number to \a n.
   */
  void assign(const int n) {v = (long double)n; }

  /**
   * Multiply the number with \a n.
   */
  void multiply(const int n) {v *= (long double)n; }

  /**
   * Print the number in the file stream \a fp.
   */
  size_t print(FILE* const fp) const {return fprintf(fp, "%Lg", v); }

  int tostring(char **str) const {
    int size=static_cast<int>( (std::log(std::abs(v))/std::log(10.0))+4 );
    *str=igraph_Calloc(size, char );
    if (! *str) {
      IGRAPH_ERROR("Cannot convert big number to string", IGRAPH_ENOMEM);
    }
    std::stringstream ss;
    ss << v;
    strncpy(*str, ss.str().c_str(), size);
    return 0;
  }
};

#endif

} //namespace bliss

#endif
