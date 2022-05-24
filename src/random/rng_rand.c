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

typedef struct {
    unsigned long int x;
} igraph_i_rng_rand_state_t;

static igraph_uint_t igraph_rng_rand_get(void *vstate) {
    igraph_i_rng_rand_state_t *state = vstate;
    state->x = (1103515245 * state->x + 12345) & 0x7fffffffUL;
    return state->x;
}

static igraph_real_t igraph_rng_rand_get_real(void *vstate) {
    return igraph_rng_rand_get (vstate) / 2147483648.0 ;
}

static igraph_error_t igraph_rng_rand_seed(void *vstate, igraph_uint_t seed) {
    igraph_i_rng_rand_state_t *state = vstate;
    state->x = (unsigned long) seed;
    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_rng_rand_init(void **state) {
    igraph_i_rng_rand_state_t *st;

    st = IGRAPH_CALLOC(1, igraph_i_rng_rand_state_t);
    if (!st) {
        IGRAPH_ERROR("Cannot initialize BSD rand RNG", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    (*state) = st;

    igraph_rng_rand_seed(st, 0);

    return IGRAPH_SUCCESS;
}

static void igraph_rng_rand_destroy(void *vstate) {
    igraph_i_rng_rand_state_t *state =
        (igraph_i_rng_rand_state_t*) vstate;
    IGRAPH_FREE(state);
}

/**
 * \var igraph_rngtype_rand
 * \brief The old BSD rand/srand random number generator.
 *
 * The sequence is
 *     <code>x_{n+1} = (a x_n + c) mod m</code>
 * with <code>a = 1103515245</code>, <code>c = 12345</code> and
 * <code>m = 2^31 = 2147483648</code>.
 * The seed specifies the initial value, <code>x_1</code>.
 *
 * </para><para>
 * The theoretical value of <code>x_{10001}</code> is 1910041713.
 *
 * </para><para>
 * The period of this generator is 2^31.
 *
 * </para><para>
 * This generator is not very good—the low bits of successive
 * numbers are correlated.
 *
 * </para><para>
 * This generator was ported from the GNU Scientific Library.
 */

const igraph_rng_type_t igraph_rngtype_rand = {
    /* name= */      "RAND",
    /* min=  */      0,
    /* max=  */      0x7fffffffUL,
    /* init= */      igraph_rng_rand_init,
    /* destroy= */   igraph_rng_rand_destroy,
    /* seed= */      igraph_rng_rand_seed,
    /* get= */       igraph_rng_rand_get,
    /* get_real= */  igraph_rng_rand_get_real,
    /* get_norm= */  0,
    /* get_geom= */  0,
    /* get_binom= */ 0,
    /* get_exp= */   0,
    /* get_gamma= */ 0
};
