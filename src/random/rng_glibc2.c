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
    int i, j;
    long int x[31];
} igraph_i_rng_glibc2_state_t;

static unsigned long int igraph_i_rng_glibc2_get(int *i, int *j, int n, long int *x) {
    unsigned long int k;

    /* The original implementation used x[*i] += x[*j] here. Considering that
     * x is signed, this is undefined behaviour according to the C standard.
     * Therefore, we temporarily cast to unsigned long int to achieve what the
     * original intention was */
    x[*i] = ((unsigned long int)x[*i]) + ((unsigned long int)x[*j]);
    k = (x[*i] >> 1) & 0x7FFFFFFF;

    (*i)++;
    if (*i == n) {
        *i = 0;
    }

    (*j)++ ;
    if (*j == n) {
        *j = 0;
    }

    return k;
}

static igraph_uint_t igraph_rng_glibc2_get(void *vstate) {
    igraph_i_rng_glibc2_state_t *state =
        (igraph_i_rng_glibc2_state_t*) vstate;
    return igraph_i_rng_glibc2_get(&state->i, &state->j, 31, state->x);
}

/* this function is independent of the bit size */

static void igraph_i_rng_glibc2_init(long int *x, int n,
                                     unsigned long int s) {
    int i;

    if (s == 0) {
        s = 1;
    }

    x[0] = (long) s;
    for (i = 1 ; i < n ; i++) {
        const long int h = s / 127773;
        const long int t = 16807 * ((long) s - h * 127773) - h * 2836;
        if (t < 0) {
            s = (unsigned long) t + 2147483647 ;
        } else {
            s = (unsigned long) t ;
        }

        x[i] = s ;
    }
}

static igraph_error_t igraph_rng_glibc2_seed(void *vstate, igraph_uint_t seed) {
    igraph_i_rng_glibc2_state_t *state =
        (igraph_i_rng_glibc2_state_t*) vstate;
    int i;

    igraph_i_rng_glibc2_init(state->x, 31, (unsigned long) seed);

    state->i = 3;
    state->j = 0;

    for (i = 0; i < 10 * 31; i++) {
        igraph_rng_glibc2_get(state);
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_rng_glibc2_init(void **state) {
    igraph_i_rng_glibc2_state_t *st;

    st = IGRAPH_CALLOC(1, igraph_i_rng_glibc2_state_t);
    IGRAPH_CHECK_OOM(st, "Cannot initialize GNU libc 2 RNG.");
    (*state) = st;

    igraph_rng_glibc2_seed(st, 0);

    return IGRAPH_SUCCESS;
}

static void igraph_rng_glibc2_destroy(void *vstate) {
    igraph_i_rng_glibc2_state_t *state =
        (igraph_i_rng_glibc2_state_t*) vstate;
    IGRAPH_FREE(state);
}

/**
 * \var igraph_rngtype_glibc2
 * \brief The random number generator introduced in GNU libc 2.
 *
 * This is a linear feedback shift register generator with a 128-byte
 * buffer. This generator was the default prior to igraph version 0.6,
 * at least on systems relying on GNU libc.
 *
 * This generator was ported from the GNU Scientific Library. It is a
 * reimplementation and does not call the system glibc generator.
 */

const igraph_rng_type_t igraph_rngtype_glibc2 = {
    /* name= */      "LIBC",
    /* bits=  */     31,
    /* init= */      igraph_rng_glibc2_init,
    /* destroy= */   igraph_rng_glibc2_destroy,
    /* seed= */      igraph_rng_glibc2_seed,
    /* get= */       igraph_rng_glibc2_get,
    /* get_int= */   0,
    /* get_real= */  0,
    /* get_norm= */  0,
    /* get_geom= */  0,
    /* get_binom= */ 0,
    /* get_exp= */   0,
    /* get_gamma= */ 0,
    /* get_pois= */  0
};
