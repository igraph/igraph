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

#include "config.h"

/* The original implementation of the 64-bit PCG random number generator in this
 * file was obtained from https://github.com/imneme/pcg-c
 *
 * PCG is dual-licensed under Apache-2.0 and MIT Licenses. MIT is compatible
 * with igraph's GPLv2 license. License notices for PCG are to be found in the
 * pcg_variants.h header
 */

#if IGRAPH_INTEGER_SIZE == 64 && defined(HAVE___UINT128_T)

#include "pcg/pcg_variants.h"

static const pcg64_random_t pcg64_initializer = PCG64_INITIALIZER;

static igraph_uint_t igraph_rng_pcg64_get(void *vstate) {
    pcg64_random_t *state = (pcg64_random_t*) vstate;
    return pcg64_random_r(state);
}

static igraph_error_t igraph_rng_pcg64_seed(void *vstate, igraph_uint_t seed) {
    pcg64_random_t *state = (pcg64_random_t*) vstate;

    if (seed == 0) {
        seed = (pcg64_initializer.inc >> 1);
    }

    /* PCG64 is seeded by a 128-bit state and a 128-bit sequence number (well, only
     * 63 bits are used from the sequence number, though). Since the unified
     * igraph RNG seeding interface provides a single igraph_uint_t as the seed,
     * we use the seed to fill in the sequence number and use the state from
     * PCG64_INITIALIZER */
    pcg64_srandom_r(state, pcg64_initializer.state, seed);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_rng_pcg64_init(void **state) {
    pcg64_random_t *st;

    st = IGRAPH_CALLOC(1, pcg64_random_t);
    IGRAPH_CHECK_OOM(st, "Cannot initialize PCG64 RNG.");
    (*state) = st;

    igraph_rng_pcg64_seed(st, 0);

    return IGRAPH_SUCCESS;
}

static void igraph_rng_pcg64_destroy(void *vstate) {
    pcg64_random_t *state = (pcg64_random_t*) vstate;
    IGRAPH_FREE(state);
}

#else

/* Dummy implementation if the compiler does not support __uint128_t */

static igraph_uint_t igraph_rng_pcg64_get(void *vstate) {
    IGRAPH_UNUSED(vstate);
    return 0;
}

static igraph_error_t igraph_rng_pcg64_seed(void *vstate, igraph_uint_t seed) {
    IGRAPH_UNUSED(vstate); IGRAPH_UNUSED(seed);
    IGRAPH_ERROR("64-bit PCG generator needs __uint128_t.", IGRAPH_UNIMPLEMENTED);
}

static igraph_error_t igraph_rng_pcg64_init(void **state) {
    IGRAPH_UNUSED(state);
    IGRAPH_ERROR("64-bit PCG generator needs __uint128_t.", IGRAPH_UNIMPLEMENTED);
}

static void igraph_rng_pcg64_destroy(void *vstate) {
    IGRAPH_UNUSED(vstate);
}

#endif

/**
 * \var igraph_rngtype_pcg64
 * \brief The PCG random number generator (64-bit version).
 *
 * This is an implementation of the PCG random number generator; see
 * https://www.pcg-random.org for more details. This implementation returns
 * 64 random bits in a single iteration. It is only available on 64-bit plaforms
 * with compilers that provide the __uint128_t type.
 *
 * </para><para>
 * PCG64 typically provides better performance than PCG32 when sampling floating
 * point numbers or very large integers, as it can provide twice as many random
 * bits in a single generation round.
 *
 * </para><para>
 * The generator was ported from the original source code published by the
 * authors at https://github.com/imneme/pcg-c.
 */

const igraph_rng_type_t igraph_rngtype_pcg64 = {
    /* name= */      "PCG64",
    /* bits=  */     64,
    /* init= */      igraph_rng_pcg64_init,
    /* destroy= */   igraph_rng_pcg64_destroy,
    /* seed= */      igraph_rng_pcg64_seed,
    /* get= */       igraph_rng_pcg64_get,
    /* get_int= */   0,
    /* get_real= */  0,
    /* get_norm= */  0,
    /* get_geom= */  0,
    /* get_binom= */ 0,
    /* get_exp= */   0,
    /* get_gamma= */ 0,
    /* get_pois= */  0
};
