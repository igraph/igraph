/* round2n.c (round floating-point number to nearest power of two) */

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
*  round2n - round floating-point number to nearest power of two
*
*  SYNOPSIS
*
*  #include "misc.h"
*  double round2n(double x);
*
*  RETURNS
*
*  Given a positive floating-point value x the routine round2n returns
*  2^n such that |x - 2^n| is minimal.
*
*  EXAMPLES
*
*  round2n(10.1) = 2^3 = 8
*  round2n(15.3) = 2^4 = 16
*  round2n(0.01) = 2^(-7) = 0.0078125
*
*  BACKGROUND
*
*  Let x = f * 2^e, where 0.5 <= f < 1 is a normalized fractional part,
*  e is an integer exponent. Then, obviously, 0.5 * 2^e <= x < 2^e, so
*  if x - 0.5 * 2^e <= 2^e - x, we choose 0.5 * 2^e = 2^(e-1), and 2^e
*  otherwise. The latter condition can be written as 2 * x <= 1.5 * 2^e
*  or 2 * f * 2^e <= 1.5 * 2^e or, finally, f <= 0.75. */

double round2n(double x)
{     int e;
      double f;
      xassert(x > 0.0);
      f = frexp(x, &e);
      return ldexp(1.0, f <= 0.75 ? e-1 : e);
}

/* eof */
