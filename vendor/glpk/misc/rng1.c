/* rng1.c (pseudo-random number generator) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2003 Free Software Foundation, Inc.
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
#include "rng.h"

/***********************************************************************
*  NAME
*
*  rng_unif_01 - obtain pseudo-random number in the range [0, 1]
*
*  SYNOPSIS
*
*  #include "rng.h"
*  double rng_unif_01(RNG *rand);
*
*  RETURNS
*
*  The routine rng_unif_01 returns a next pseudo-random number which is
*  uniformly distributed in the range [0, 1]. */

double rng_unif_01(RNG *rand)
{     double x;
      x = (double)rng_next_rand(rand) / 2147483647.0;
      xassert(0.0 <= x && x <= 1.0);
      return x;
}

/***********************************************************************
*  NAME
*
*  rng_uniform - obtain pseudo-random number in the range [a, b]
*
*  SYNOPSIS
*
*  #include "rng.h"
*  double rng_uniform(RNG *rand, double a, double b);
*
*  RETURNS
*
*  The routine rng_uniform returns a next pseudo-random number which is
*  uniformly distributed in the range [a, b]. */

double rng_uniform(RNG *rand, double a, double b)
{     double x;
      xassert(a < b);
      x = rng_unif_01(rand);
      x = a * (1.0 - x) + b * x;
      xassert(a <= x && x <= b);
      return x;
}

/* eof */
