/* fp2rat.c (convert floating-point number to rational number) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2000-2013 Free Software Foundation, Inc.
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

#include "env.h"
#include "misc.h"

/***********************************************************************
*  NAME
*
*  fp2rat - convert floating-point number to rational number
*
*  SYNOPSIS
*
*  #include "misc.h"
*  int fp2rat(double x, double eps, double *p, double *q);
*
*  DESCRIPTION
*
*  Given a floating-point number 0 <= x < 1 the routine fp2rat finds
*  its "best" rational approximation p / q, where p >= 0 and q > 0 are
*  integer numbers, such that |x - p / q| <= eps.
*
*  RETURNS
*
*  The routine fp2rat returns the number of iterations used to achieve
*  the specified precision eps.
*
*  EXAMPLES
*
*  For x = sqrt(2) - 1 = 0.414213562373095 and eps = 1e-6 the routine
*  gives p = 408 and q = 985, where 408 / 985 = 0.414213197969543.
*
*  BACKGROUND
*
*  It is well known that every positive real number x can be expressed
*  as the following continued fraction:
*
*     x = b[0] + a[1]
*                ------------------------
*                b[1] + a[2]
*                       -----------------
*                       b[2] + a[3]
*                              ----------
*                              b[3] + ...
*
*  where:
*
*     a[k] = 1,                  k = 0, 1, 2, ...
*
*     b[k] = floor(x[k]),        k = 0, 1, 2, ...
*
*     x[0] = x,
*
*     x[k] = 1 / frac(x[k-1]),   k = 1, 2, 3, ...
*
*  To find the "best" rational approximation of x the routine computes
*  partial fractions f[k] by dropping after k terms as follows:
*
*     f[k] = A[k] / B[k],
*
*  where:
*
*     A[-1] = 1,   A[0] = b[0],   B[-1] = 0,   B[0] = 1,
*
*     A[k] = b[k] * A[k-1] + a[k] * A[k-2],
*
*     B[k] = b[k] * B[k-1] + a[k] * B[k-2].
*
*  Once the condition
*
*     |x - f[k]| <= eps
*
*  has been satisfied, the routine reports p = A[k] and q = B[k] as the
*  final answer.
*
*  In the table below here is some statistics obtained for one million
*  random numbers uniformly distributed in the range [0, 1).
*
*      eps      max p   mean p      max q    mean q  max k   mean k
*     -------------------------------------------------------------
*     1e-1          8      1.6          9       3.2    3      1.4
*     1e-2         98      6.2         99      12.4    5      2.4
*     1e-3        997     20.7        998      41.5    8      3.4
*     1e-4       9959     66.6       9960     133.5   10      4.4
*     1e-5      97403    211.7      97404     424.2   13      5.3
*     1e-6     479669    669.9     479670    1342.9   15      6.3
*     1e-7    1579030   2127.3    3962146    4257.8   16      7.3
*     1e-8   26188823   6749.4   26188824   13503.4   19      8.2
*
*  REFERENCES
*
*  W. B. Jones and W. J. Thron, "Continued Fractions: Analytic Theory
*  and Applications," Encyclopedia on Mathematics and Its Applications,
*  Addison-Wesley, 1980. */

int fp2rat(double x, double eps, double *p, double *q)
{     int k;
      double xk, Akm1, Ak, Bkm1, Bk, ak, bk, fk, temp;
      xassert(0.0 <= x && x < 1.0);
      for (k = 0; ; k++)
      {  xassert(k <= 100);
         if (k == 0)
         {  /* x[0] = x */
            xk = x;
            /* A[-1] = 1 */
            Akm1 = 1.0;
            /* A[0] = b[0] = floor(x[0]) = 0 */
            Ak = 0.0;
            /* B[-1] = 0 */
            Bkm1 = 0.0;
            /* B[0] = 1 */
            Bk = 1.0;
         }
         else
         {  /* x[k] = 1 / frac(x[k-1]) */
            temp = xk - floor(xk);
            xassert(temp != 0.0);
            xk = 1.0 / temp;
            /* a[k] = 1 */
            ak = 1.0;
            /* b[k] = floor(x[k]) */
            bk = floor(xk);
            /* A[k] = b[k] * A[k-1] + a[k] * A[k-2] */
            temp = bk * Ak + ak * Akm1;
            Akm1 = Ak, Ak = temp;
            /* B[k] = b[k] * B[k-1] + a[k] * B[k-2] */
            temp = bk * Bk + ak * Bkm1;
            Bkm1 = Bk, Bk = temp;
         }
         /* f[k] = A[k] / B[k] */
         fk = Ak / Bk;
#if 0
         print("%.*g / %.*g = %.*g",
            DBL_DIG, Ak, DBL_DIG, Bk, DBL_DIG, fk);
#endif
         if (fabs(x - fk) <= eps)
            break;
      }
      *p = Ak;
      *q = Bk;
      return k;
}

/* eof */
