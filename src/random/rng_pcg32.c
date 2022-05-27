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

/* The original implementation of the 32-bit PCG random number generator in this
 * file was downloaded from https://www.pcg-random.org/download.html
 *
 * PCG is licensed under the Apache License. The original license notice follows
 * below.
 *
 * Some modifications were made to the generator to make it fit igraph's RNG
 * interface. Functions related to the global PCG RNG and bounded RNG were
 * removed as we don't need a global RNG.
 */

/*
 * PCG Random Number Generation for C.
 *
 * Copyright 2014 Melissa O'Neill <oneill@pcg-random.org>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * For additional information about the PCG random number generation scheme,
 * including its license and other licensing options, visit
 *
 *       http://www.pcg-random.org
 */

struct pcg_state_setseq_64 {    // Internals are *Private*.
    uint64_t state;             // RNG state.  All values are possible.
    uint64_t inc;               // Controls which RNG sequence (stream) is
                                // selected. Must *always* be odd.
};
typedef struct pcg_state_setseq_64 pcg32_random_t;

// If you *must* statically initialize it, here's one.

#define PCG32_INITIALIZER   { 0x853c49e6748fea9bULL, 0xda3e39cb94b95bdbULL }

static uint32_t pcg32_random_r(pcg32_random_t* rng);

// pcg32_srandom_r(rng, initstate, initseq):
//     Seed the rng.  Specified in two parts, state initializer and a
//     sequence selection constant (a.k.a. stream id)

static void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate, uint64_t initseq)
{
    rng->state = 0U;
    rng->inc = (initseq << 1u) | 1u;
    pcg32_random_r(rng);
    rng->state += initstate;
    pcg32_random_r(rng);
}

// pcg32_random_r(rng)
//     Generate a uniformly distributed 32-bit random number

static uint32_t pcg32_random_r(pcg32_random_t* rng)
{
    uint64_t oldstate = rng->state;
    rng->state = oldstate * 6364136223846793005ULL + rng->inc;
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

/* igraph-specific code follows below */

static const pcg32_random_t pcg32_initializer = PCG32_INITIALIZER;

static igraph_uint_t igraph_rng_pcg32_get(void *vstate) {
    pcg32_random_t *state = (pcg32_random_t*) vstate;
    return pcg32_random_r(state);
}

static igraph_error_t igraph_rng_pcg32_seed(void *vstate, igraph_uint_t seed) {
    pcg32_random_t *state = (pcg32_random_t*) vstate;

    if (seed == 0) {
        seed = (pcg32_initializer.inc >> 1);
    }

    /* PCG is seeded by a 64-bit state and a 64-bit sequence number (well, only
     * 63 bits are used from the sequence number, though). Since the unified
     * igraph RNG seeding interface provides a single igraph_uint_t as the seed,
     * we use the seed to fill in the sequence number and use the state from
     * PCG32_INITIALIZER */
    pcg32_srandom_r(state, pcg32_initializer.state, seed);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_rng_pcg32_init(void **state) {
    pcg32_random_t *st;

    st = IGRAPH_CALLOC(1, pcg32_random_t);
    IGRAPH_CHECK_OOM(st, "Cannot initialize PCG32 RNG");
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
 * \brief The PCG random number generator.
 *
 * This is an implementation of the PCG random number generator; see
 * https://www.pcg-random.org for more details. This implementation returns
 * 32 random bits in a single iteration.
 *
 * The generator was ported from the original source code published by the
 * authors at https://www.pcg-random.org .
 */

const igraph_rng_type_t igraph_rngtype_pcg32 = {
    /* name= */      "PCG32",
    /* bits=  */     32,
    /* init= */      igraph_rng_pcg32_init,
    /* destroy= */   igraph_rng_pcg32_destroy,
    /* seed= */      igraph_rng_pcg32_seed,
    /* get= */       igraph_rng_pcg32_get,
    /* get_real= */  0,
    /* get_norm= */  0,
    /* get_geom= */  0,
    /* get_binom= */ 0,
    /* get_exp= */   0,
    /* get_gamma= */ 0
};
