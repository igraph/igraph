/* gcd.c (greatest common divisor) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2000 Free Software Foundation, Inc.
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
*  gcd - find greatest common divisor of two integers
*
*  SYNOPSIS
*
*  #include "misc.h"
*  int gcd(int x, int y);
*
*  RETURNS
*
*  The routine gcd returns gcd(x, y), the greatest common divisor of
*  the two positive integers given.
*
*  ALGORITHM
*
*  The routine gcd is based on Euclid's algorithm.
*
*  REFERENCES
*
*  Don Knuth, The Art of Computer Programming, Vol.2: Seminumerical
*  Algorithms, 3rd Edition, Addison-Wesley, 1997. Section 4.5.2: The
*  Greatest Common Divisor, pp. 333-56. */

int gcd(int x, int y)
{     int r;
      xassert(x > 0 && y > 0);
      while (y > 0)
         r = x % y, x = y, y = r;
      return x;
}

/***********************************************************************
*  NAME
*
*  gcdn - find greatest common divisor of n integers
*
*  SYNOPSIS
*
*  #include "misc.h"
*  int gcdn(int n, int x[]);
*
*  RETURNS
*
*  The routine gcdn returns gcd(x[1], x[2], ..., x[n]), the greatest
*  common divisor of n positive integers given, n > 0.
*
*  BACKGROUND
*
*  The routine gcdn is based on the following identity:
*
*     gcd(x, y, z) = gcd(gcd(x, y), z).
*
*  REFERENCES
*
*  Don Knuth, The Art of Computer Programming, Vol.2: Seminumerical
*  Algorithms, 3rd Edition, Addison-Wesley, 1997. Section 4.5.2: The
*  Greatest Common Divisor, pp. 333-56. */

int gcdn(int n, int x[])
{     int d, j;
      xassert(n > 0);
      for (j = 1; j <= n; j++)
      {  xassert(x[j] > 0);
         if (j == 1)
            d = x[1];
         else
            d = gcd(d, x[j]);
         if (d == 1)
            break;
      }
      return d;
}

/* eof */
