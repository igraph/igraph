/* glprng.h (pseudo-random number generator) */

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

#ifndef GLPRNG_H
#define GLPRNG_H

typedef struct RNG RNG;

struct RNG
{     /* Knuth's portable pseudo-random number generator */
      int A[56];
      /* pseudo-random values */
      int *fptr;
      /* the next A value to be exported */
};

#define rng_create_rand _glp_rng_create_rand
RNG *rng_create_rand(void);
/* create pseudo-random number generator */

#define rng_init_rand _glp_rng_init_rand
void rng_init_rand(RNG *rand, int seed);
/* initialize pseudo-random number generator */

#define rng_next_rand _glp_rng_next_rand
int rng_next_rand(RNG *rand);
/* obtain pseudo-random integer in the range [0, 2^31-1] */

#define rng_unif_rand _glp_rng_unif_rand
int rng_unif_rand(RNG *rand, int m);
/* obtain pseudo-random integer in the range [0, m-1] */

#define rng_delete_rand _glp_rng_delete_rand
void rng_delete_rand(RNG *rand);
/* delete pseudo-random number generator */

#define rng_unif_01 _glp_rng_unif_01
double rng_unif_01(RNG *rand);
/* obtain pseudo-random number in the range [0, 1] */

#define rng_uniform _glp_rng_uniform
double rng_uniform(RNG *rand, double a, double b);
/* obtain pseudo-random number in the range [a, b] */

#endif

/* eof */
