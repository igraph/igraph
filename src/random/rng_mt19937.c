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

#include <string.h> /* memset() */
#include <ctype.h>


#define N 624   /* Period parameters */
#define M 397

/* most significant w-r bits */
static const uint32_t UPPER_MASK = UINT32_C(0x80000000);

/* least significant r bits */
static const uint32_t LOWER_MASK = UINT32_C(0x7fffffff);

typedef struct {
    uint32_t mt[N];
    int mti;
} igraph_i_rng_mt19937_state_t;

static igraph_uint_t igraph_rng_mt19937_get(void *vstate) {
    igraph_i_rng_mt19937_state_t *state = vstate;

    uint32_t k;
    uint32_t *const mt = state->mt;

#define MAGIC(y) (((y) & 0x1) ? UINT32_C(0x9908b0df) : 0)

    if (state->mti >= N) {
        /* generate N words at one time */
        int kk;

        for (kk = 0; kk < N - M; kk++) {
            uint32_t y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + M] ^ (y >> 1) ^ MAGIC(y);
        }
        for (; kk < N - 1; kk++) {
            uint32_t y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ MAGIC(y);
        }

        {
            uint32_t y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
            mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ MAGIC(y);
        }

        state->mti = 0;
    }

#undef MAGIC

    /* Tempering */

    k = mt[state->mti];
    k ^= (k >> 11);
    k ^= (k << 7) & UINT32_C(0x9d2c5680);
    k ^= (k << 15) & UINT32_C(0xefc60000);
    k ^= (k >> 18);

    state->mti++;

    return k;
}

static igraph_error_t igraph_rng_mt19937_seed(void *vstate, igraph_uint_t seed) {
    igraph_i_rng_mt19937_state_t *state = vstate;
    int i;

    memset(state, 0, sizeof(igraph_i_rng_mt19937_state_t));

    if (seed == 0) {
        seed = 4357;   /* the default seed is 4357 */
    }
    state->mt[0] = seed & UINT32_C(0xffffffff);

    for (i = 1; i < N; i++) {
        /* See Knuth's "Art of Computer Programming" Vol. 2, 3rd
           Ed. p.106 for multiplier. */
        state->mt[i] =
            (UINT32_C(1812433253) * (state->mt[i - 1] ^ (state->mt[i - 1] >> 30)) +
             (uint32_t) i);
        state->mt[i] &= UINT32_C(0xffffffff);
    }

    state->mti = i;
    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_rng_mt19937_init(void **state) {
    igraph_i_rng_mt19937_state_t *st;

    st = IGRAPH_CALLOC(1, igraph_i_rng_mt19937_state_t);
    IGRAPH_CHECK_OOM(st, "Cannot initialize MT19937 RNG.");
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
 * \c igraph_rng_mt19937_seed() reproduces this. Later versions switched to
 * 5489 as the default seed, you can choose this explicitly via
 * \ref igraph_rng_seed() instead if you require it.
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
    /* get_int= */   0,
    /* get_real= */  0,
    /* get_norm= */  0,
    /* get_geom= */  0,
    /* get_binom= */ 0,
    /* get_exp= */   0,
    /* get_gamma= */ 0,
    /* get_pois= */  0
};

#undef N
#undef M
