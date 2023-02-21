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

#include "pcg/pcg_variants.h"

#include "config.h" /* IGRAPH_THREAD_LOCAL */

/* The original implementation of the 32-bit PCG random number generator in this
 * file was obtained from https://github.com/imneme/pcg-c
 *
 * PCG is dual-licensed under Apache-2.0 and MIT Licenses. MIT is compatible
 * with igraph's GPLv2 license. License notices for PCG are to be found in the
 * pcg_variants.h header
 */

static const pcg32_random_t pcg32_initializer = PCG32_INITIALIZER;

static igraph_uint_t igraph_rng_pcg32_get(void *vstate) {
    pcg32_random_t *state = (pcg32_random_t*) vstate;
    return pcg32_random_r(state);
}

static igraph_error_t igraph_rng_pcg32_seed(void *vstate, igraph_uint_t seed) {
    pcg32_random_t *state = (pcg32_random_t*) vstate;

    /* PCG32 is seeded by a 64-bit state and a 64-bit sequence number (well, only
     * 63 bits are used from the sequence number, though). Since the unified
     * igraph RNG seeding interface provides a single igraph_uint_t as the seed,
     * we use the seed to fill in the sequence number and use the state from
     * PCG32_INITIALIZER */
    if (seed == 0) {
        /* If you feel the temptation to unify the two branches by running
         * seed = pcg32_initializer.inc >> 1, don't.
         * seed is an igraph_uint_t, so it can be 32-bit or 64-bit.
         * pcg32_initializer.inc is always 64-bit.
         */
        pcg32_srandom_r(state, pcg32_initializer.state, pcg32_initializer.inc >> 1);
    } else {
        pcg32_srandom_r(state, pcg32_initializer.state, seed);
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_rng_pcg32_init(void **state) {
    pcg32_random_t *st;

    st = IGRAPH_CALLOC(1, pcg32_random_t);
    IGRAPH_CHECK_OOM(st, "Cannot initialize PCG32 RNG.");
    (*state) = st;

    igraph_rng_pcg32_seed(st, 0);

    return IGRAPH_SUCCESS;
}

static void igraph_rng_pcg32_destroy(void *vstate) {
    pcg32_random_t *state = (pcg32_random_t*) vstate;
    IGRAPH_FREE(state);
}

/**
 * \var igraph_rngtype_pcg32
 * \brief The PCG random number generator (32-bit version).
 *
 * This is an implementation of the PCG random number generator; see
 * https://www.pcg-random.org for more details. This implementation returns
 * 32 random bits in a single iteration.
 *
 * </para><para>
 * The generator was ported from the original source code published by the
 * authors at https://github.com/imneme/pcg-c.
 */

const igraph_rng_type_t igraph_rngtype_pcg32 = {
    /* name= */      "PCG32",
    /* bits=  */     32,
    /* init= */      igraph_rng_pcg32_init,
    /* destroy= */   igraph_rng_pcg32_destroy,
    /* seed= */      igraph_rng_pcg32_seed,
    /* get= */       igraph_rng_pcg32_get,
    /* get_int= */   0,
    /* get_real= */  0,
    /* get_norm= */  0,
    /* get_geom= */  0,
    /* get_binom= */ 0,
    /* get_exp= */   0,
    /* get_gamma= */ 0,
    /* get_pois= */  0
};

/***** Default RNG, used upon igraph startup *****/

#define addr(a) (&a)

static pcg32_random_t igraph_i_rng_default_state = PCG32_INITIALIZER;

IGRAPH_THREAD_LOCAL igraph_rng_t igraph_i_rng_default = {
    addr(igraph_rngtype_pcg32),
    addr(igraph_i_rng_default_state),
    /* is_seeded = */ true
};

#undef addr
