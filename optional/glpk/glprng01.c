/* glprng01.c */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  This code is a modified version of the module GB_FLIP, a portable
*  pseudo-random number generator. The original version of GB_FLIP is
*  a part of The Stanford GraphBase developed by Donald E. Knuth (see
*  http://www-cs-staff.stanford.edu/~knuth/sgb.html).
*
*  Note that all changes concern only external names, so this modified
*  version produces exactly the same results as the original version.
*
*  Changes were made by Andrew Makhorin <mao@gnu.org>.
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
#include "glprng.h"

#if 0
int A[56] = { -1 };
#else
#define A (rand->A)
#endif
/* pseudo-random values */

#if 0
int *fptr = A;
#else
#define fptr (rand->fptr)
#endif
/* the next A value to be exported */

#define mod_diff(x, y) (((x) - (y)) & 0x7FFFFFFF)
/* difference modulo 2^31 */

static int flip_cycle(RNG *rand)
{     /* this is an auxiliary routine to do 55 more steps of the basic
         recurrence, at high speed, and to reset fptr */
      int *ii, *jj;
      for (ii = &A[1], jj = &A[32]; jj <= &A[55]; ii++, jj++)
         *ii = mod_diff(*ii, *jj);
      for (jj = &A[1]; ii <= &A[55]; ii++, jj++)
         *ii = mod_diff(*ii, *jj);
      fptr = &A[54];
      return A[55];
}

/***********************************************************************
*  NAME
*
*  rng_create_rand - create pseudo-random number generator
*
*  SYNOPSIS
*
*  #include "glprng.h"
*  RNG *rng_create_rand(void);
*
*  DESCRIPTION
*
*  The routine rng_create_rand creates and initializes a pseudo-random
*  number generator.
*
*  RETURNS
*
*  The routine returns a pointer to the generator created. */

RNG *rng_create_rand(void)
{     RNG *rand;
      int i;
      rand = xmalloc(sizeof(RNG));
      A[0] = -1;
      for (i = 1; i <= 55; i++) A[i] = 0;
      fptr = A;
      rng_init_rand(rand, 1);
      return rand;
}

/***********************************************************************
*  NAME
*
*  rng_init_rand - initialize pseudo-random number generator
*
*  SYNOPSIS
*
*  #include "glprng.h"
*  void rng_init_rand(RNG *rand, int seed);
*
*  DESCRIPTION
*
*  The routine rng_init_rand initializes the pseudo-random number
*  generator. The parameter seed may be any integer number. Note that
*  on creating the generator this routine is called with the parameter
*  seed equal to 1. */

void rng_init_rand(RNG *rand, int seed)
{     int i;
      int prev = seed, next = 1;
      seed = prev = mod_diff(prev, 0);
      A[55] = prev;
      for (i = 21; i; i = (i + 21) % 55)
      {  A[i] = next;
         next = mod_diff(prev, next);
         if (seed & 1)
            seed = 0x40000000 + (seed >> 1);
         else
            seed >>= 1;
         next = mod_diff(next, seed);
         prev = A[i];
      }
      flip_cycle(rand);
      flip_cycle(rand);
      flip_cycle(rand);
      flip_cycle(rand);
      flip_cycle(rand);
      return;
}

/***********************************************************************
*  NAME
*
*  rng_next_rand - obtain pseudo-random integer in the range [0, 2^31-1]
*
*  SYNOPSIS
*
*  #include "glprng.h"
*  int rng_next_rand(RNG *rand);
*
*  RETURNS
*
*  The routine rng_next_rand returns a next pseudo-random integer which
*  is uniformly distributed between 0 and 2^31-1, inclusive. The period
*  length of the generated numbers is 2^85 - 2^30. The low order bits of
*  the generated numbers are just as random as the high-order bits. */

int rng_next_rand(RNG *rand)
{     return
         *fptr >= 0 ? *fptr-- : flip_cycle(rand);
}

/***********************************************************************
*  NAME
*
*  rng_unif_rand - obtain pseudo-random integer in the range [0, m-1]
*
*  SYNOPSIS
*
*  #include "glprng.h"
*  int rng_unif_rand(RNG *rand, int m);
*
*  RETURNS
*
*  The routine rng_unif_rand returns a next pseudo-random integer which
*  is uniformly distributed between 0 and m-1, inclusive, where m is any
*  positive integer less than 2^31. */

#define two_to_the_31 ((unsigned int)0x80000000)

int rng_unif_rand(RNG *rand, int m)
{     unsigned int t = two_to_the_31 - (two_to_the_31 % m);
      int r;
      xassert(m > 0);
      do { r = rng_next_rand(rand); } while (t <= (unsigned int)r);
      return r % m;
}

/***********************************************************************
*  NAME
*
*  rng_delete_rand - delete pseudo-random number generator
*
*  SYNOPSIS
*
*  #include "glprng.h"
*  void rng_delete_rand(RNG *rand);
*
*  DESCRIPTION
*
*  The routine rng_delete_rand frees all the memory allocated to the
*  specified pseudo-random number generator. */

void rng_delete_rand(RNG *rand)
{     xfree(rand);
      return;
}

/**********************************************************************/

#if 0
/* To be sure that this modified version produces the same results as
   the original version, run this validation program. */

int main(void)
{     RNG *rand;
      int j;
      rand = rng_create_rand();
      rng_init_rand(rand, -314159);
      if (rng_next_rand(rand) != 119318998)
      {  fprintf(stderr, "Failure on the first try!\n");
         return -1;
      }
      for (j = 1; j <= 133; j++) rng_next_rand(rand);
      if (rng_unif_rand(rand, 0x55555555) != 748103812)
      {  fprintf(stderr, "Failure on the second try!\n");
         return -2;
      }
      fprintf(stderr, "OK, the random-number generator routines seem to"
         " work!\n");
      rng_delete_rand(rand);
      return 0;
}
#endif

/* eof */
