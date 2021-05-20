/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

#include "test_utilities.inc"

/* How many "true" elements in the Boolean vector? */
long vector_bool_count(const igraph_vector_bool_t *vec) {
    long i, n = igraph_vector_bool_size(vec), cnt = 0;

    for (i=0; i < n; ++i) {
        if (VECTOR(*vec)[i]) {
            cnt++;
        }
    }
    return cnt;
}

int main() {
    igraph_t g;
    igraph_vector_t type_dist;
    igraph_matrix_t pref_mat;
    igraph_vector_t types, in_types, out_types;
    igraph_bool_t connected, has_loop, has_multi;
    igraph_vector_bool_t is_loop;
    long i;

    igraph_vector_init(&types, 0);

    /* Symmetric preference game */
    igraph_vector_init_real(&type_dist, 3, 1.0, 1.0, 1.0);

    igraph_matrix_init(&pref_mat, 3, 3);
    for (i = 0; i < 3; i++) {
        MATRIX(pref_mat, i, i) = 0.2;
    }

    /* undirected, no loops */
    IGRAPH_CHECK(igraph_preference_game(&g, 1000, 3, &type_dist, /*fixed_sizes=*/ 0,
                                        &pref_mat, &types, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS));

    IGRAPH_ASSERT(igraph_vcount(&g) == 1000);
    IGRAPH_ASSERT(! igraph_is_directed(&g));

    igraph_is_connected(&g, &connected, IGRAPH_STRONG);
    IGRAPH_ASSERT(! connected);

    igraph_has_loop(&g, &has_loop);
    IGRAPH_ASSERT(! has_loop);

    igraph_has_multiple(&g, &has_multi);
    IGRAPH_ASSERT(! has_multi);

    IGRAPH_ASSERT(igraph_vector_size(&types) == igraph_vcount(&g));
    IGRAPH_ASSERT(igraph_vector_min(&types) == 0);
    IGRAPH_ASSERT(igraph_vector_max(&types) == 2);

    igraph_destroy(&g);

    /* Note: preference matrix must be symmetric in the undirected case. */
    for (i = 0; i < 2; i++) {
        MATRIX(pref_mat, i, i + 1) = 0.1;
        MATRIX(pref_mat, i + 1, i) = 0.1;
    }

    /* directed, no loops */
    IGRAPH_CHECK(igraph_preference_game(&g, 1000, 3, &type_dist, /*fixed_sizes=*/0,
                                        &pref_mat, &types, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS));

    IGRAPH_ASSERT(igraph_vcount(&g) == 1000);
    IGRAPH_ASSERT(igraph_is_directed(&g));

    igraph_has_loop(&g, &has_loop);
    IGRAPH_ASSERT(! has_loop);

    igraph_has_multiple(&g, &has_multi);
    IGRAPH_ASSERT(! has_multi);

    IGRAPH_ASSERT(igraph_vector_size(&types) == igraph_vcount(&g));
    IGRAPH_ASSERT(igraph_vector_min(&types) == 0);
    IGRAPH_ASSERT(igraph_vector_max(&types) == 2);

    igraph_destroy(&g);

    /* undirected, loops */
    for (i = 0; i < 3; i++) {
        MATRIX(pref_mat, i, i) = 1.0;
    }

    IGRAPH_CHECK(igraph_preference_game(&g, 100, 3, &type_dist, /*fixed_sizes=*/ 0,
                                        &pref_mat, &types, IGRAPH_UNDIRECTED, IGRAPH_LOOPS));

    IGRAPH_ASSERT(igraph_vcount(&g) == 100);
    IGRAPH_ASSERT(igraph_ecount(&g) >= 1395);
    IGRAPH_ASSERT(!igraph_is_directed(&g));

    igraph_has_loop(&g, &has_loop);
    IGRAPH_ASSERT(has_loop);

    igraph_has_multiple(&g, &has_multi);
    IGRAPH_ASSERT(! has_multi);

    IGRAPH_ASSERT(igraph_vector_size(&types) == igraph_vcount(&g));
    IGRAPH_ASSERT(igraph_vector_min(&types) == 0);
    IGRAPH_ASSERT(igraph_vector_max(&types) == 2);

    igraph_destroy(&g);

    /* directed, loops */
    IGRAPH_CHECK(igraph_preference_game(&g, 100, 3, &type_dist, /*fixed_sizes=*/ 0,
                                        &pref_mat, NULL, IGRAPH_DIRECTED, IGRAPH_LOOPS));

    IGRAPH_ASSERT(igraph_vcount(&g) == 100);
    IGRAPH_ASSERT(igraph_ecount(&g) >= 2700);
    IGRAPH_ASSERT(igraph_is_directed(&g));

    igraph_has_loop(&g, &has_loop);
    IGRAPH_ASSERT(has_loop);

    igraph_has_multiple(&g, &has_multi);
    IGRAPH_ASSERT(! has_multi);

    igraph_destroy(&g);

    igraph_vector_destroy(&types);

    /* Asymmetric preference game */

    igraph_vector_init(&in_types, 0);
    igraph_vector_init(&out_types, 0);

    /* directed, no loops */
    igraph_matrix_resize(&pref_mat, 2, 3);
    MATRIX(pref_mat, 0, 0) = 1;
    MATRIX(pref_mat, 0, 1) = 1;
    MATRIX(pref_mat, 0, 2) = 1;
    MATRIX(pref_mat, 1, 0) = 1;
    MATRIX(pref_mat, 1, 1) = 1;
    MATRIX(pref_mat, 1, 2) = 1;

    IGRAPH_CHECK(igraph_asymmetric_preference_game(&g, 100, 2, 3, NULL, &pref_mat, &in_types, &out_types, IGRAPH_NO_LOOPS));

    IGRAPH_ASSERT(igraph_vcount(&g) == 100);
    IGRAPH_ASSERT(igraph_ecount(&g) == 9900);
    IGRAPH_ASSERT(igraph_is_directed(&g));

    igraph_has_loop(&g, &has_loop);
    IGRAPH_ASSERT(! has_loop);

    igraph_has_multiple(&g, &has_multi);
    IGRAPH_ASSERT(! has_multi);

    igraph_destroy(&g);

    /* directed, loops */
    igraph_matrix_resize(&pref_mat, 2, 2);
    MATRIX(pref_mat, 0, 0) = 1;
    MATRIX(pref_mat, 0, 1) = 1;
    MATRIX(pref_mat, 1, 0) = 1;
    MATRIX(pref_mat, 1, 1) = 1;

    IGRAPH_CHECK(igraph_asymmetric_preference_game(&g, 100, 2, 2, NULL, &pref_mat, NULL, NULL, IGRAPH_LOOPS));

    IGRAPH_ASSERT(igraph_vcount(&g) == 100);
    IGRAPH_ASSERT(igraph_ecount(&g) == 10000);
    IGRAPH_ASSERT(igraph_is_directed(&g));

    igraph_vector_bool_init(&is_loop, 0);
    igraph_is_loop(&g, &is_loop, igraph_ess_all(IGRAPH_EDGEORDER_ID));
    IGRAPH_ASSERT(vector_bool_count(&is_loop) == 100);
    igraph_vector_bool_destroy(&is_loop);

    igraph_has_multiple(&g, &has_multi);
    IGRAPH_ASSERT(! has_multi);

    igraph_destroy(&g);

    igraph_vector_destroy(&type_dist);
    igraph_matrix_destroy(&pref_mat);

    igraph_vector_destroy(&in_types);
    igraph_vector_destroy(&out_types);

    VERIFY_FINALLY_STACK();

    return 0;
}
