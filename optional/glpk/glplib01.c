/* glplib01.c (bignum arithmetic) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008,
*  2009, 2010 Andrew Makhorin, Department for Applied Informatics,
*  Moscow Aviation Institute, Moscow, Russia. All rights reserved.
*  E-mail: <mao@gnu.org>.
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

#include "glpenv.h"
#include "glplib.h"

/***********************************************************************
*  Two routines below are intended to multiply and divide unsigned
*  integer numbers of arbitrary precision.
* 
*  The routines assume that an unsigned integer number is represented in
*  the positional numeral system with the base 2^16 = 65536, i.e. each
*  "digit" of the number is in the range [0, 65535] and represented as
*  a 16-bit value of the unsigned short type. In other words, a number x
*  has the following representation:
* 
*         n-1
*     x = sum d[j] * 65536^j,
*         j=0
* 
*  where n is the number of places (positions), and d[j] is j-th "digit"
*  of x, 0 <= d[j] <= 65535.
***********************************************************************/

/***********************************************************************
*  NAME
*
*  bigmul - multiply unsigned integer numbers of arbitrary precision
* 
*  SYNOPSIS
* 
*  #include "glplib.h"
*  void bigmul(int n, int m, unsigned short x[], unsigned short y[]);
* 
*  DESCRIPTION
* 
*  The routine bigmul multiplies unsigned integer numbers of arbitrary
*  precision.
* 
*  n is the number of digits of multiplicand, n >= 1;
* 
*  m is the number of digits of multiplier, m >= 1;
* 
*  x is an array containing digits of the multiplicand in elements
*  x[m], x[m+1], ..., x[n+m-1]. Contents of x[0], x[1], ..., x[m-1] are
*  ignored on entry.
* 
*  y is an array containing digits of the multiplier in elements y[0],
*  y[1], ..., y[m-1].
* 
*  On exit digits of the product are stored in elements x[0], x[1], ...,
*  x[n+m-1]. The array y is not changed. */

void bigmul(int n, int m, unsigned short x[], unsigned short y[])
{     int i, j;
      unsigned int t;
      xassert(n >= 1);
      xassert(m >= 1);
      for (j = 0; j < m; j++) x[j] = 0;
      for (i = 0; i < n; i++)
      {  if (x[i+m])
         {  t = 0;
            for (j = 0; j < m; j++)
            {  t += (unsigned int)x[i+m] * (unsigned int)y[j] +
                    (unsigned int)x[i+j];
               x[i+j] = (unsigned short)t;
               t >>= 16;
            }
            x[i+m] = (unsigned short)t;
         }
      }
      return;
}

/***********************************************************************
*  NAME
*
*  bigdiv - divide unsigned integer numbers of arbitrary precision
* 
*  SYNOPSIS
* 
*  #include "glplib.h"
*  void bigdiv(int n, int m, unsigned short x[], unsigned short y[]);
* 
*  DESCRIPTION
* 
*  The routine bigdiv divides one unsigned integer number of arbitrary
*  precision by another with the algorithm described in [1].
* 
*  n is the difference between the number of digits of dividend and the
*  number of digits of divisor, n >= 0.
* 
*  m is the number of digits of divisor, m >= 1.
* 
*  x is an array containing digits of the dividend in elements x[0],
*  x[1], ..., x[n+m-1].
* 
*  y is an array containing digits of the divisor in elements y[0],
*  y[1], ..., y[m-1]. The highest digit y[m-1] must be non-zero.
* 
*  On exit n+1 digits of the quotient are stored in elements x[m],
*  x[m+1], ..., x[n+m], and m digits of the remainder are stored in
*  elements x[0], x[1], ..., x[m-1]. The array y is changed but then
*  restored.
*
*  REFERENCES
*
*  1. D. Knuth. The Art of Computer Programming. Vol. 2: Seminumerical
*  Algorithms. Stanford University, 1969. */

void bigdiv(int n, int m, unsigned short x[], unsigned short y[])
{     int i, j;
      unsigned int t;
      unsigned short d, q, r;
      xassert(n >= 0);
      xassert(m >= 1);
      xassert(y[m-1] != 0);
      /* special case when divisor has the only digit */
      if (m == 1)
      {  d = 0;
         for (i = n; i >= 0; i--)
         {  t = ((unsigned int)d << 16) + (unsigned int)x[i];
            x[i+1] = (unsigned short)(t / y[0]);
            d = (unsigned short)(t % y[0]);
         }
         x[0] = d;
         goto done;
      }
      /* multiply dividend and divisor by a normalizing coefficient in
         order to provide the condition y[m-1] >= base / 2 */
      d = (unsigned short)(0x10000 / ((unsigned int)y[m-1] + 1));
      if (d == 1)
         x[n+m] = 0;
      else
      {  t = 0;
         for (i = 0; i < n+m; i++)
         {  t += (unsigned int)x[i] * (unsigned int)d;
            x[i] = (unsigned short)t;
            t >>= 16;
         }
         x[n+m] = (unsigned short)t;
         t = 0;
         for (j = 0; j < m; j++)
         {  t += (unsigned int)y[j] * (unsigned int)d;
            y[j] = (unsigned short)t;
            t >>= 16;
         }
      }
      /* main loop */
      for (i = n; i >= 0; i--)
      {  /* estimate and correct the current digit of quotient */
         if (x[i+m] < y[m-1])
         {  t = ((unsigned int)x[i+m] << 16) + (unsigned int)x[i+m-1];
            q = (unsigned short)(t / (unsigned int)y[m-1]);
            r = (unsigned short)(t % (unsigned int)y[m-1]);
            if (q == 0) goto putq; else goto test;
         }
         q = 0;
         r = x[i+m-1];
decr:    q--; /* if q = 0 then q-- = 0xFFFF */
         t = (unsigned int)r + (unsigned int)y[m-1];
         r = (unsigned short)t;
         if (t > 0xFFFF) goto msub;
test:    t = (unsigned int)y[m-2] * (unsigned int)q;
         if ((unsigned short)(t >> 16) > r) goto decr;
         if ((unsigned short)(t >> 16) < r) goto msub;
         if ((unsigned short)t > x[i+m-2]) goto decr;
msub:    /* now subtract divisor multiplied by the current digit of
            quotient from the current dividend */
         if (q == 0) goto putq;
         t = 0;
         for (j = 0; j < m; j++)
         {  t += (unsigned int)y[j] * (unsigned int)q;
            if (x[i+j] < (unsigned short)t) t += 0x10000;
            x[i+j] -= (unsigned short)t;
            t >>= 16;
         }
         if (x[i+m] >= (unsigned short)t) goto putq;
         /* perform correcting addition, because the current digit of
            quotient is greater by one than its correct value */
         q--;
         t = 0;
         for (j = 0; j < m; j++)
         {  t += (unsigned int)x[i+j] + (unsigned int)y[j];
            x[i+j] = (unsigned short)t;
            t >>= 16;
         }
putq:    /* store the current digit of quotient */
         x[i+m] = q;
      }
      /* divide divisor and remainder by the normalizing coefficient in
         order to restore their original values */
      if (d > 1)
      {  t = 0;
         for (i = m-1; i >= 0; i--)
         {  t = (t << 16) + (unsigned int)x[i];
            x[i] = (unsigned short)(t / (unsigned int)d);
            t %= (unsigned int)d;
         }
         t = 0;
         for (j = m-1; j >= 0; j--)
         {  t = (t << 16) + (unsigned int)y[j];
            y[j] = (unsigned short)(t / (unsigned int)d);
            t %= (unsigned int)d;
         }
      }
done: return;
}

/**********************************************************************/

#if 0
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "glprng.h"

#define N_MAX 7
/* maximal number of digits in multiplicand */

#define M_MAX 5
/* maximal number of digits in multiplier */

#define N_TEST 1000000
/* number of tests */

int main(void)
{     RNG *rand;
      int d, j, n, m, test;
      unsigned short x[N_MAX], y[M_MAX], z[N_MAX+M_MAX];
      rand = rng_create_rand();
      for (test = 1; test <= N_TEST; test++)
      {  /* x[0,...,n-1] := multiplicand */
         n = 1 + rng_unif_rand(rand, N_MAX-1);
         assert(1 <= n && n <= N_MAX);
         for (j = 0; j < n; j++)
         {  d = rng_unif_rand(rand, 65536);
            assert(0 <= d && d <= 65535);
            x[j] = (unsigned short)d;
         }
         /* y[0,...,m-1] := multiplier */
         m = 1 + rng_unif_rand(rand, M_MAX-1);
         assert(1 <= m && m <= M_MAX);
         for (j = 0; j < m; j++)
         {  d = rng_unif_rand(rand, 65536);
            assert(0 <= d && d <= 65535);
            y[j] = (unsigned short)d;
         }
         if (y[m-1] == 0) y[m-1] = 1;
         /* z[0,...,n+m-1] := x * y */
         for (j = 0; j < n; j++) z[m+j] = x[j];
         bigmul(n, m, z, y);
         /* z[0,...,m-1] := z mod y, z[m,...,n+m-1] := z div y */
         bigdiv(n, m, z, y);
         /* z mod y must be 0 */
         for (j = 0; j < m; j++) assert(z[j] == 0);
         /* z div y must be x */
         for (j = 0; j < n; j++) assert(z[m+j] == x[j]);
      }
      fprintf(stderr, "%d tests successfully passed\n", N_TEST);
      rng_delete_rand(rand);
      return 0;
}
#endif

/* eof */
