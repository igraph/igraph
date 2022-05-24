/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "igraph_random.h"

#include "igraph_memory.h"
#include "igraph_types.h"

#include "config.h" /* IGRAPH_THREAD_LOCAL */

#include <string.h> /* memset() */


#define N 624   /* Period parameters */
#define M 397

/* most significant w-r bits */
static const unsigned long UPPER_MASK = 0x80000000UL;

/* least significant r bits */
static const unsigned long LOWER_MASK = 0x7fffffffUL;

typedef struct {
    unsigned long mt[N];
    int mti;
} igraph_i_rng_mt19937_state_t;

static igraph_uint_t igraph_rng_mt19937_get(void *vstate) {
    igraph_i_rng_mt19937_state_t *state = vstate;

    unsigned long k;
    unsigned long int *const mt = state->mt;

#define MAGIC(y) (((y)&0x1) ? 0x9908b0dfUL : 0)

    if (state->mti >= N) {
        /* generate N words at one time */
        int kk;

        for (kk = 0; kk < N - M; kk++) {
            unsigned long y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + M] ^ (y >> 1) ^ MAGIC(y);
        }
        for (; kk < N - 1; kk++) {
            unsigned long y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ MAGIC(y);
        }

        {
            unsigned long y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
            mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ MAGIC(y);
        }

        state->mti = 0;
    }

#undef MAGIC

    /* Tempering */

    k = mt[state->mti];
    k ^= (k >> 11);
    k ^= (k << 7) & 0x9d2c5680UL;
    k ^= (k << 15) & 0xefc60000UL;
    k ^= (k >> 18);

    state->mti++;

    return k;
}

static igraph_real_t igraph_rng_mt19937_get_real(void *vstate) {
    return igraph_rng_mt19937_get (vstate) / 4294967296.0 ;
}

static igraph_error_t igraph_rng_mt19937_seed(void *vstate, igraph_uint_t seed) {
    igraph_i_rng_mt19937_state_t *state = vstate;
    int i;

    memset(state, 0, sizeof(igraph_i_rng_mt19937_state_t));

    if (seed == 0) {
        seed = 4357;   /* the default seed is 4357 */
    }
    state->mt[0] = seed & 0xffffffffUL;

    for (i = 1; i < N; i++) {
        /* See Knuth's "Art of Computer Programming" Vol. 2, 3rd
           Ed. p.106 for multiplier. */
        state->mt[i] =
            (1812433253UL * (state->mt[i - 1] ^ (state->mt[i - 1] >> 30)) +
             (unsigned long) i);
        state->mt[i] &= 0xffffffffUL;
    }

    state->mti = i;
    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_rng_mt19937_init(void **state) {
    igraph_i_rng_mt19937_state_t *st;

    st = IGRAPH_CALLOC(1, igraph_i_rng_mt19937_state_t);
    if (!st) {
        IGRAPH_ERROR("Cannot initialize MT19937 RNG", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    (*state) = st;

    igraph_rng_mt19937_seed(st, 0);

    return IGRAPH_SUCCESS;
}

static void igraph_rng_mt19937_destroy(void *vstate) {
    igraph_i_rng_mt19937_state_t *state =
        (igraph_i_rng_mt19937_state_t*) vstate;
    IGRAPH_FREE(state);
}

/**
 * \var igraph_rngtype_mt19937
 * \brief The MT19937 random number generator.
 *
 * The MT19937 generator of Makoto Matsumoto and Takuji Nishimura is a
 * variant of the twisted generalized feedback shift-register
 * algorithm, and is known as the “Mersenne Twister” generator. It has
 * a Mersenne prime period of 2^19937 - 1 (about 10^6000) and is
 * equi-distributed in 623 dimensions. It has passed the diehard
 * statistical tests. It uses 624 words of state per generator and is
 * comparable in speed to the other generators. The original generator
 * used a default seed of 4357 and choosing \c s equal to zero in
 * \c gsl_rng_set reproduces this. Later versions switched to 5489 as the
 * default seed, you can choose this explicitly via \ref igraph_rng_seed()
 * instead if you require it.
 *
 * </para><para>
 * For more information see,
 * Makoto Matsumoto and Takuji Nishimura, “Mersenne Twister: A
 * 623-dimensionally equidistributed uniform pseudorandom number
 * generator”. ACM Transactions on Modeling and Computer Simulation,
 * Vol. 8, No. 1 (Jan. 1998), Pages 3–30
 *
 * </para><para>
 * The generator \c igraph_rngtype_mt19937 uses the second revision of the
 * seeding procedure published by the two authors above in 2002. The
 * original seeding procedures could cause spurious artifacts for some
 * seed values.
 *
 * </para><para>
 * This generator was ported from the GNU Scientific Library.
 */

const igraph_rng_type_t igraph_rngtype_mt19937 = {
    /* name= */      "MT19937",
    /* bits=  */     32,
    /* init= */      igraph_rng_mt19937_init,
    /* destroy= */   igraph_rng_mt19937_destroy,
    /* seed= */      igraph_rng_mt19937_seed,
    /* get= */       igraph_rng_mt19937_get,
    /* get_real= */  igraph_rng_mt19937_get_real,
    /* get_norm= */  0,
    /* get_geom= */  0,
    /* get_binom= */ 0,
    /* get_exp= */   0,
    /* get_gamma= */ 0
};

#undef N
#undef M


/***** Default RNG, used upon igraph startup *****/

#define addr(a) (&a)

static igraph_i_rng_mt19937_state_t igraph_i_rng_default_state;

IGRAPH_THREAD_LOCAL igraph_rng_t igraph_i_rng_default = {
    addr(igraph_rngtype_mt19937),
    addr(igraph_i_rng_default_state),
    /* def= */ 1
};

#undef addr
