/* glplib02.c (64-bit arithmetic) */

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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wsometimes-uninitialized"
#endif

#include "glpenv.h"
#include "glplib.h"

/***********************************************************************
*  NAME
*
*  xlset - expand integer to long integer
*
*  SYNOPSIS
*
*  #include "glplib.h"
*  glp_long xlset(int x);
*
*  RETURNS
*
*  The routine xlset returns x expanded to long integer. */

glp_long xlset(int x)
{     glp_long t;
      t.lo = x, t.hi = (x >= 0 ? 0 : -1);
      return t;
}

/***********************************************************************
*  NAME
*
*  xlneg - negate long integer
*
*  SYNOPSIS
*
*  #include "glplib.h"
*  glp_long xlneg(glp_long x);
*
*  RETURNS
*
*  The routine xlneg returns the difference  0 - x. */

glp_long xlneg(glp_long x)
{     if (x.lo)
         x.lo = - x.lo, x.hi = ~x.hi;
      else
         x.hi = - x.hi;
      return x;
}

/***********************************************************************
*  NAME
*
*  xladd - add long integers
*
*  SYNOPSIS
*
*  #include "glplib.h"
*  glp_long xladd(glp_long x, glp_long y);
*
*  RETURNS
*
*  The routine xladd returns the sum x + y. */

glp_long xladd(glp_long x, glp_long y)
{     if ((unsigned int)x.lo <= 0xFFFFFFFF - (unsigned int)y.lo)
         x.lo += y.lo, x.hi += y.hi;
      else
         x.lo += y.lo, x.hi += y.hi + 1;
      return x;
}

/***********************************************************************
*  NAME
*
*  xlsub - subtract long integers
*
*  SYNOPSIS
*
*  #include "glplib.h"
*  glp_long xlsub(glp_long x, glp_long y);
*
*  RETURNS
*
*  The routine xlsub returns the difference x - y. */

glp_long xlsub(glp_long x, glp_long y)
{     return
         xladd(x, xlneg(y));
}

/***********************************************************************
*  NAME
*
*  xlcmp - compare long integers
*
*  SYNOPSIS
*
*  #include "glplib.h"
*  int xlcmp(glp_long x, glp_long y);
*
*  RETURNS
*
*  The routine xlcmp returns the sign of the difference x - y. */

int xlcmp(glp_long x, glp_long y)
{     if (x.hi >= 0 && y.hi <  0) return +1;
      if (x.hi <  0 && y.hi >= 0) return -1;
      if ((unsigned int)x.hi < (unsigned int)y.hi) return -1;
      if ((unsigned int)x.hi > (unsigned int)y.hi) return +1;
      if ((unsigned int)x.lo < (unsigned int)y.lo) return -1;
      if ((unsigned int)x.lo > (unsigned int)y.lo) return +1;
      return 0;
}

/***********************************************************************
*  NAME
*
*  xlmul - multiply long integers
*
*  SYNOPSIS
*
*  #include "glplib.h"
*  glp_long xlmul(glp_long x, glp_long y);
*
*  RETURNS
*
*  The routine xlmul returns the product x * y. */

glp_long xlmul(glp_long x, glp_long y)
{     unsigned short xx[8], yy[4];
      xx[4] = (unsigned short)x.lo;
      xx[5] = (unsigned short)(x.lo >> 16);
      xx[6] = (unsigned short)x.hi;
      xx[7] = (unsigned short)(x.hi >> 16);
      yy[0] = (unsigned short)y.lo;
      yy[1] = (unsigned short)(y.lo >> 16);
      yy[2] = (unsigned short)y.hi;
      yy[3] = (unsigned short)(y.hi >> 16);
      bigmul(4, 4, xx, yy);
      x.lo = (unsigned int)xx[0] | ((unsigned int)xx[1] << 16);
      x.hi = (unsigned int)xx[2] | ((unsigned int)xx[3] << 16);
      return x;
}

/***********************************************************************
*  NAME
*
*  xldiv - divide long integers
*
*  SYNOPSIS
*
*  #include "glplib.h"
*  glp_ldiv xldiv(glp_long x, glp_long y);
*
*  RETURNS
*
*  The routine xldiv returns a structure of type glp_ldiv containing
*  members quot (the quotient) and rem (the remainder), both of type
*  glp_long. */

glp_ldiv xldiv(glp_long x, glp_long y)
{     glp_ldiv t;
      int m, sx, sy;
      unsigned short xx[8], yy[4];
      /* sx := sign(x) */
      sx = (x.hi < 0);
      /* sy := sign(y) */
      sy = (y.hi < 0);
      /* x := |x| */
      if (sx) x = xlneg(x);
      /* y := |y| */
      if (sy) y = xlneg(y);
      /* compute x div y and x mod y */
      xx[0] = (unsigned short)x.lo;
      xx[1] = (unsigned short)(x.lo >> 16);
      xx[2] = (unsigned short)x.hi;
      xx[3] = (unsigned short)(x.hi >> 16);
      yy[0] = (unsigned short)y.lo;
      yy[1] = (unsigned short)(y.lo >> 16);
      yy[2] = (unsigned short)y.hi;
      yy[3] = (unsigned short)(y.hi >> 16);
      if (yy[3])
         m = 4;
      else if (yy[2])
         m = 3;
      else if (yy[1])
         m = 2;
      else if (yy[0])
         m = 1;
      else
         xerror("xldiv: divide by zero\n");
      bigdiv(4 - m, m, xx, yy);
      /* remainder in x[0], x[1], ..., x[m-1] */
      t.rem.lo = (unsigned int)xx[0], t.rem.hi = 0;
      if (m >= 2) t.rem.lo |= (unsigned int)xx[1] << 16;
      if (m >= 3) t.rem.hi = (unsigned int)xx[2];
      if (m >= 4) t.rem.hi |= (unsigned int)xx[3] << 16;
      if (sx) t.rem = xlneg(t.rem);
      /* quotient in x[m], x[m+1], ..., x[4] */
      t.quot.lo = (unsigned int)xx[m], t.quot.hi = 0;
      if (m <= 3) t.quot.lo |= (unsigned int)xx[m+1] << 16;
      if (m <= 2) t.quot.hi = (unsigned int)xx[m+2];
      if (m <= 1) t.quot.hi |= (unsigned int)xx[m+3] << 16;
      if (sx ^ sy) t.quot = xlneg(t.quot);
      return t;
}

/***********************************************************************
*  NAME
*
*  xltod - convert long integer to double
*
*  SYNOPSIS
*
*  #include "glplib.h"
*  double xltod(glp_long x);
*
*  RETURNS
*
*  The routine xltod returns x converted to double. */

double xltod(glp_long x)
{     double s, z;
      if (x.hi >= 0)
         s = +1.0;
      else
         s = -1.0, x = xlneg(x);
      if (x.hi >= 0)
         z = 4294967296.0 * (double)x.hi + (double)(unsigned int)x.lo;
      else
      {  xassert(x.hi == 0x80000000 && x.lo == 0x00000000);
         z = 9223372036854775808.0; /* 2^63 */
      }
      return s * z;
}

char *xltoa(glp_long x, char *s)
{     /* convert long integer to character string */
      static const char *d = "0123456789";
      glp_ldiv t;
      int neg, len;
      if (x.hi >= 0)
         neg = 0;
      else
         neg = 1, x = xlneg(x);
      if (x.hi >= 0)
      {  len = 0;
         while (!(x.hi == 0 && x.lo == 0))
         {  t = xldiv(x, xlset(10));
            xassert(0 <= t.rem.lo && t.rem.lo <= 9);
            s[len++] = d[t.rem.lo];
            x = t.quot;
         }
         if (len == 0) s[len++] = d[0];
         if (neg) s[len++] = '-';
         s[len] = '\0';
         strrev(s);
      }
      else
         strcpy(s, "-9223372036854775808"); /* -2^63 */
      return s;
}

/**********************************************************************/

#if 0
#include "glprng.h"

#define N_TEST 1000000
/* number of tests */

static glp_long myrand(RNG *rand)
{     glp_long x;
      int k;
      k = rng_unif_rand(rand, 4);
      xassert(0 <= k && k <= 3);
      x.lo = rng_unif_rand(rand, 65536);
      if (k == 1 || k == 3)
      {  x.lo <<= 16;
         x.lo += rng_unif_rand(rand, 65536);
      }
      if (k <= 1)
         x.hi = 0;
      else
         x.hi = rng_unif_rand(rand, 65536);
      if (k == 3)
      {  x.hi <<= 16;
         x.hi += rng_unif_rand(rand, 65536);
      }
      if (rng_unif_rand(rand, 2)) x = xlneg(x);
      return x;
}

int main(void)
{     RNG *rand;
      glp_long x, y;
      glp_ldiv z;
      int test;
      rand = rng_create_rand();
      for (test = 1; test <= N_TEST; test++)
      {  x = myrand(rand);
         y = myrand(rand);
         if (y.lo == 0 && y.hi == 0) y.lo = 1;
         /* z.quot := x div y, z.rem := x mod y */
         z = xldiv(x, y);
         /* x must be equal to y * z.quot + z.rem */
         xassert(xlcmp(x, xladd(xlmul(y, z.quot), z.rem)) == 0);
      }
      xprintf("%d tests successfully passed\n", N_TEST);
      rng_delete_rand(rand);
      return 0;
}
#endif

/* eof */
