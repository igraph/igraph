/* mygmp.h (integer and rational arithmetic) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2008-2015 Free Software Foundation, Inc.
*  Written by Andrew Makhorin <mao@gnu.org>.
*
*  GLPK is free software: you can redistribute it and/or modify it
*  under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  GLPK is distributed in the hope that it will be useful, but WITHOUT
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
*  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
*  License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with GLPK. If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/

#ifndef MYGMP_H
#define MYGMP_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_GMP               /* use GNU MP library */

#include <gmp.h>

#define gmp_pool_count() 0

#define gmp_free_mem() ((void)0)

#else                         /* use GLPK MP module */

/***********************************************************************
*  INTEGER NUMBERS
*  ---------------
*  Depending on its magnitude an integer number of arbitrary precision
*  is represented either in short format or in long format.
*
*  Short format corresponds to the int type and allows representing
*  integer numbers in the range [-(2^31-1), +(2^31-1)]. Note that for
*  the most negative number of int type the short format is not used.
*
*  In long format integer numbers are represented using the positional
*  system with the base (radix) 2^16 = 65536:
*
*     x = (-1)^s sum{j in 0..n-1} d[j] * 65536^j,
*
*  where x is the integer to be represented, s is its sign (+1 or -1),
*  d[j] are its digits (0 <= d[j] <= 65535).
*
*  RATIONAL NUMBERS
*  ----------------
*  A rational number is represented as an irreducible fraction:
*
*     p / q,
*
*  where p (numerator) and q (denominator) are integer numbers (q > 0)
*  having no common divisors. */

struct mpz
{     /* integer number */
      int val;
      /* if ptr is a null pointer, the number is in short format, and
         val is its value; otherwise, the number is in long format, and
         val is its sign (+1 or -1) */
      struct mpz_seg *ptr;
      /* pointer to the linked list of the number segments ordered in
         ascending of powers of the base */
};

struct mpz_seg
{     /* integer number segment */
      unsigned short d[6];
      /* six digits of the number ordered in ascending of powers of the
         base */
      struct mpz_seg *next;
      /* pointer to the next number segment */
};

struct mpq
{     /* rational number (p / q) */
      struct mpz p;
      /* numerator */
      struct mpz q;
      /* denominator */
};

typedef struct mpz *mpz_t;
typedef struct mpq *mpq_t;

#define gmp_get_atom _glp_gmp_get_atom
void *gmp_get_atom(int size);

#define gmp_free_atom _glp_gmp_free_atom
void gmp_free_atom(void *ptr, int size);

#define gmp_pool_count _glp_gmp_pool_count
int gmp_pool_count(void);

#define gmp_get_work _glp_gmp_get_work
unsigned short *gmp_get_work(int size);

#define gmp_free_mem _glp_gmp_free_mem
void gmp_free_mem(void);

#define mpz_init(x) (void)((x) = _mpz_init())

#define _mpz_init _glp_mpz_init
mpz_t _mpz_init(void);
/* initialize x and set its value to 0 */

#define mpz_clear _glp_mpz_clear
void mpz_clear(mpz_t x);
/* free the space occupied by x */

#define mpz_set _glp_mpz_set
void mpz_set(mpz_t z, mpz_t x);
/* set the value of z from x */

#define mpz_set_si _glp_mpz_set_si
void mpz_set_si(mpz_t x, int val);
/* set the value of x to val */

#define mpz_get_d _glp_mpz_get_d
double mpz_get_d(mpz_t x);
/* convert x to a double, truncating if necessary */

#define mpz_get_d_2exp _glp_mpz_get_d_2exp
double mpz_get_d_2exp(int *exp, mpz_t x);
/* convert x to a double, returning the exponent separately */

#define mpz_swap _glp_mpz_swap
void mpz_swap(mpz_t x, mpz_t y);
/* swap the values x and y efficiently */

#define mpz_add _glp_mpz_add
void mpz_add(mpz_t, mpz_t, mpz_t);
/* set z to x + y */

#define mpz_sub _glp_mpz_sub
void mpz_sub(mpz_t, mpz_t, mpz_t);
/* set z to x - y */

#define mpz_mul _glp_mpz_mul
void mpz_mul(mpz_t, mpz_t, mpz_t);
/* set z to x * y */

#define mpz_neg _glp_mpz_neg
void mpz_neg(mpz_t z, mpz_t x);
/* set z to 0 - x */

#define mpz_abs _glp_mpz_abs
void mpz_abs(mpz_t z, mpz_t x);
/* set z to the absolute value of x */

#define mpz_div _glp_mpz_div
void mpz_div(mpz_t q, mpz_t r, mpz_t x, mpz_t y);
/* divide x by y, forming quotient q and/or remainder r */

#define mpz_gcd _glp_mpz_gcd
void mpz_gcd(mpz_t z, mpz_t x, mpz_t y);
/* set z to the greatest common divisor of x and y */

#define mpz_cmp _glp_mpz_cmp
int mpz_cmp(mpz_t x, mpz_t y);
/* compare x and y */

#define mpz_sgn _glp_mpz_sgn
int mpz_sgn(mpz_t x);
/* return +1 if x > 0, 0 if x = 0, and -1 if x < 0 */

#define mpz_out_str _glp_mpz_out_str
int mpz_out_str(void *fp, int base, mpz_t x);
/* output x on stream fp, as a string in given base */

#define mpq_init(x) (void)((x) = _mpq_init())

#define _mpq_init _glp_mpq_init
mpq_t _mpq_init(void);
/* initialize x, and set its value to 0/1 */

#define mpq_clear _glp_mpq_clear
void mpq_clear(mpq_t x);
/* free the space occupied by x */

#define mpq_canonicalize _glp_mpq_canonicalize
void mpq_canonicalize(mpq_t x);
/* canonicalize x */

#define mpq_set _glp_mpq_set
void mpq_set(mpq_t z, mpq_t x);
/* set the value of z from x */

#define mpq_set_si _glp_mpq_set_si
void mpq_set_si(mpq_t x, int p, unsigned int q);
/* set the value of x to p/q */

#define mpq_get_d _glp_mpq_get_d
double mpq_get_d(mpq_t x);
/* convert x to a double, truncating if necessary */

#define mpq_set_d _glp_mpq_set_d
void mpq_set_d(mpq_t x, double val);
/* set x to val; there is no rounding, the conversion is exact */

#define mpq_add _glp_mpq_add
void mpq_add(mpq_t z, mpq_t x, mpq_t y);
/* set z to x + y */

#define mpq_sub _glp_mpq_sub
void mpq_sub(mpq_t z, mpq_t x, mpq_t y);
/* set z to x - y */

#define mpq_mul _glp_mpq_mul
void mpq_mul(mpq_t z, mpq_t x, mpq_t y);
/* set z to x * y */

#define mpq_div _glp_mpq_div
void mpq_div(mpq_t z, mpq_t x, mpq_t y);
/* set z to x / y */

#define mpq_neg _glp_mpq_neg
void mpq_neg(mpq_t z, mpq_t x);
/* set z to 0 - x */

#define mpq_abs _glp_mpq_abs
void mpq_abs(mpq_t z, mpq_t x);
/* set z to the absolute value of x */

#define mpq_cmp _glp_mpq_cmp
int mpq_cmp(mpq_t x, mpq_t y);
/* compare x and y */

#define mpq_sgn _glp_mpq_sgn
int mpq_sgn(mpq_t x);
/* return +1 if x > 0, 0 if x = 0, and -1 if x < 0 */

#define mpq_out_str _glp_mpq_out_str
int mpq_out_str(void *fp, int base, mpq_t x);
/* output x on stream fp, as a string in given base */

#endif

#endif

/* eof */
