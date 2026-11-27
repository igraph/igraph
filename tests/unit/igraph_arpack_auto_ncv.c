/*
   igraph library.
   Copyright (C) 2026  The igraph development team <igraph@igraph.org>

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

/* Regression test for the automatic choice of ncv in igraph_arpack_rssolve().
 * With ncv < n, ARPACK failed with "No shifts could be applied" or "No
 * eigenvalues to sufficient accuracy" on small graphs whose largest
 * eigenvalue is repeated, depending on the random seed. */

#include <igraph.h>

#include "test_utilities.h"

/* Multiplies by diag(1, 2, ..., n). */
static igraph_error_t diag_matvec(igraph_real_t *to, const igraph_real_t *from, int n, void *extra) {
    IGRAPH_UNUSED(extra);
    for (igraph_int_t i = 0; i < n; i++) {
        to[i] = (i + 1) * from[i];
    }
    return IGRAPH_SUCCESS;
}

int main(void) {
    igraph_t g, k4;
    igraph_vector_t hub, ec;
    igraph_vector_int_t membership;
    igraph_real_t value;

    /* The functions below warn about disconnected or undirected input. */
    igraph_set_warning_handler(igraph_warning_handler_ignore);

    /* Two disjoint copies of K4. The leading eigenvector method failed for
     * 14% to 19% of seeds, with both the vendored ARPACK and arpack-ng. */
    igraph_full(&k4, 4, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_disjoint_union(&g, &k4, &k4);
    igraph_destroy(&k4);

    igraph_vector_int_init(&membership, 0);
    for (igraph_int_t seed = 0; seed < 100; seed++) {
        igraph_rng_seed(igraph_rng_default(), seed);
        igraph_community_leading_eigenvector(&g, NULL, NULL, &membership, 8, NULL, NULL,
                                             false, NULL, NULL, NULL, NULL, NULL);
        /* The communities are the two cliques. */
        for (igraph_int_t i = 0; i < 8; i++) {
            IGRAPH_ASSERT(VECTOR(membership)[i] == VECTOR(membership)[i < 4 ? 0 : 4]);
        }
        IGRAPH_ASSERT(VECTOR(membership)[0] != VECTOR(membership)[4]);
    }
    igraph_vector_int_destroy(&membership);
    igraph_destroy(&g);

    /* A perfect matching on 6 vertices, the originally reported case. This
     * failed for about 1% of seeds when the vendored ARPACK was compiled
     * with fused multiply-add. */
    igraph_small(&g, 6, IGRAPH_UNDIRECTED, 1,2, 0,3, 4,5, -1);

    igraph_vector_init(&hub, 0);
    igraph_vector_init(&ec, 0);
    for (igraph_int_t seed = 0; seed < 1000; seed++) {
        igraph_rng_seed(igraph_rng_default(), seed);
        igraph_hub_and_authority_scores(&g, &hub, NULL, &value, NULL, NULL);
        IGRAPH_ASSERT(fabs(value - 1) < 1e-10);

        igraph_rng_seed(igraph_rng_default(), seed);
        igraph_eigenvector_centrality(&g, &ec, &value, IGRAPH_ALL, NULL, NULL);
        IGRAPH_ASSERT(fabs(value - 1) < 1e-10);
    }
    igraph_vector_destroy(&ec);
    igraph_vector_destroy(&hub);
    igraph_destroy(&g);

    /* The automatic ncv must fit in caller-provided storage when possible. */
    {
        igraph_arpack_options_t options;
        igraph_arpack_storage_t storage;
        igraph_vector_t values;
        igraph_matrix_t cvalues;

        igraph_vector_init(&values, 0);
        igraph_matrix_init(&cvalues, 0, 0);

        igraph_arpack_options_init(&options);
        options.n = 6;
        options.nev = 1;
        options.which[0] = 'L'; options.which[1] = 'A';
        igraph_arpack_storage_init(&storage, 6, 5, 6, true);
        igraph_arpack_rssolve(diag_matvec, NULL, &options, &storage, &values, NULL);
        IGRAPH_ASSERT(fabs(VECTOR(values)[0] - 6) < 1e-10);
        igraph_arpack_storage_destroy(&storage);

        igraph_arpack_options_init(&options);
        options.n = 30;
        options.nev = 1;
        options.which[0] = 'L'; options.which[1] = 'M';
        igraph_arpack_storage_init(&storage, 30, 10, 30, false);
        igraph_arpack_rnsolve(diag_matvec, NULL, &options, &storage, &cvalues, NULL);
        IGRAPH_ASSERT(fabs(MATRIX(cvalues, 0, 0) - 30) < 1e-10);
        igraph_arpack_storage_destroy(&storage);

        igraph_matrix_destroy(&cvalues);
        igraph_vector_destroy(&values);
    }

    VERIFY_FINALLY_STACK();

    return 0;
}
