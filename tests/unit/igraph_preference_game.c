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

#include "test_utilities.h"

/* How many "true" elements in the Boolean vector? */
igraph_integer_t vector_bool_count(const igraph_vector_bool_t *vec) {
    igraph_integer_t i, n = igraph_vector_bool_size(vec), cnt = 0;

    for (i=0; i < n; ++i) {
        if (VECTOR(*vec)[i]) {
            cnt++;
        }
    }
    return cnt;
}

int main(void) {
    igraph_t g;
    igraph_vector_t type_dist;
    igraph_matrix_t pref_mat, type_dist_mat;
    igraph_vector_int_t types, out_types, in_types;
    igraph_bool_t connected, has_loop, has_multi;
    igraph_vector_bool_t is_loop;
    igraph_integer_t i, j, count;

    igraph_vector_int_init(&types, 0);

    /* Symmetric preference game */
    igraph_vector_init_real(&type_dist, 3, 1.0, 1.0, 1.0);

    igraph_matrix_init(&pref_mat, 3, 3);
    for (i = 0; i < 3; i++) {
        MATRIX(pref_mat, i, i) = 0.2;
    }

    /* undirected, no loops */
    igraph_preference_game(&g, 1000, 3, &type_dist, /*fixed_sizes=*/ 0,
                           &pref_mat, &types, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

    IGRAPH_ASSERT(igraph_vcount(&g) == 1000);
    IGRAPH_ASSERT(! igraph_is_directed(&g));

    igraph_is_connected(&g, &connected, IGRAPH_STRONG);
    IGRAPH_ASSERT(! connected);

    igraph_has_loop(&g, &has_loop);
    IGRAPH_ASSERT(! has_loop);

    igraph_has_multiple(&g, &has_multi);
    IGRAPH_ASSERT(! has_multi);

    IGRAPH_ASSERT(igraph_vector_int_size(&types) == igraph_vcount(&g));
    IGRAPH_ASSERT(igraph_vector_int_min(&types) == 0);
    IGRAPH_ASSERT(igraph_vector_int_max(&types) == 2);

    igraph_destroy(&g);

    /* Note: preference matrix must be symmetric in the undirected case. */
    for (i = 0; i < 2; i++) {
        MATRIX(pref_mat, i, i + 1) = 0.1;
        MATRIX(pref_mat, i + 1, i) = 0.1;
    }

    /* directed, no loops */
    igraph_preference_game(&g, 1000, 3, &type_dist, /*fixed_sizes=*/0,
                           &pref_mat, &types, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);

    IGRAPH_ASSERT(igraph_vcount(&g) == 1000);
    IGRAPH_ASSERT(igraph_is_directed(&g));

    igraph_has_loop(&g, &has_loop);
    IGRAPH_ASSERT(! has_loop);

    igraph_has_multiple(&g, &has_multi);
    IGRAPH_ASSERT(! has_multi);

    IGRAPH_ASSERT(igraph_vector_int_size(&types) == igraph_vcount(&g));
    IGRAPH_ASSERT(igraph_vector_int_min(&types) == 0);
    IGRAPH_ASSERT(igraph_vector_int_max(&types) == 2);

    igraph_destroy(&g);

    /* undirected, loops */
    for (i = 0; i < 3; i++) {
        MATRIX(pref_mat, i, i) = 1.0;
    }

    igraph_preference_game(&g, 100, 3, &type_dist, /*fixed_sizes=*/ 0,
                           &pref_mat, &types, IGRAPH_UNDIRECTED, IGRAPH_LOOPS);

    IGRAPH_ASSERT(igraph_vcount(&g) == 100);
    IGRAPH_ASSERT(igraph_ecount(&g) >= 1395);
    IGRAPH_ASSERT(!igraph_is_directed(&g));

    igraph_has_loop(&g, &has_loop);
    IGRAPH_ASSERT(has_loop);

    igraph_has_multiple(&g, &has_multi);
    IGRAPH_ASSERT(! has_multi);

    IGRAPH_ASSERT(igraph_vector_int_size(&types) == igraph_vcount(&g));
    IGRAPH_ASSERT(igraph_vector_int_min(&types) == 0);
    IGRAPH_ASSERT(igraph_vector_int_max(&types) == 2);

    igraph_destroy(&g);

    /* directed, loops */
    igraph_preference_game(&g, 100, 3, &type_dist, /*fixed_sizes=*/ 0,
                           &pref_mat, NULL, IGRAPH_DIRECTED, IGRAPH_LOOPS);

    IGRAPH_ASSERT(igraph_vcount(&g) == 100);
    IGRAPH_ASSERT(igraph_ecount(&g) >= 2700);
    IGRAPH_ASSERT(igraph_is_directed(&g));

    igraph_has_loop(&g, &has_loop);
    IGRAPH_ASSERT(has_loop);

    igraph_has_multiple(&g, &has_multi);
    IGRAPH_ASSERT(! has_multi);

    igraph_destroy(&g);

    /* fixed sizes, divide evenly */
    igraph_matrix_resize(&pref_mat, 9, 9);
    for (i = 0; i < 9; i++) {
        for (j = 0; j < 9; j++) {
            MATRIX(pref_mat, i, j) = (j == i + 1 || j == i - 1) ? 0.1 : 0;
        }
    }
    igraph_preference_game(&g, 50, 9, /*type_dist=*/ 0, /*fixed_sizes=*/ 1,
                           &pref_mat, &types, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

    IGRAPH_ASSERT(igraph_vcount(&g) == 50);
    IGRAPH_ASSERT(!igraph_is_directed(&g));

    igraph_has_loop(&g, &has_loop);
    IGRAPH_ASSERT(! has_loop);

    igraph_has_multiple(&g, &has_multi);
    IGRAPH_ASSERT(! has_multi);

    for (i = 0; i < 9; i++) {
        count = 0;
        for (j = 0; j < 50; j++) {
            if (VECTOR(types)[j] == i) {
                count++;
            }
        }
        IGRAPH_ASSERT(count == 5 || count == 6);
    }

    igraph_destroy(&g);

    igraph_vector_int_destroy(&types);

    /* Asymmetric preference game */

    igraph_vector_int_init(&out_types, 0);
    igraph_vector_int_init(&in_types, 0);

    /* directed, no loops */
    igraph_matrix_resize(&pref_mat, 2, 3);
    MATRIX(pref_mat, 0, 0) = 1;
    MATRIX(pref_mat, 0, 1) = 1;
    MATRIX(pref_mat, 0, 2) = 1;
    MATRIX(pref_mat, 1, 0) = 1;
    MATRIX(pref_mat, 1, 1) = 1;
    MATRIX(pref_mat, 1, 2) = 1;

    igraph_asymmetric_preference_game(&g, 100, 2, 3, NULL, &pref_mat, &out_types, &in_types, IGRAPH_NO_LOOPS);

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

    igraph_asymmetric_preference_game(&g, 100, 2, 2, NULL, &pref_mat, NULL, NULL, IGRAPH_LOOPS);

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

    /* check that vertex out-/in-types are generated correctly; here pref_mat does not matter */

    igraph_matrix_resize(&pref_mat, 3, 2);
    igraph_matrix_null(&pref_mat);

    igraph_matrix_init(&type_dist_mat, 3, 2);
    igraph_matrix_null(&type_dist_mat);
    MATRIX(type_dist_mat, 2, 0) = 1; /* all out-types are 2, all in-types are 0 */

    igraph_asymmetric_preference_game(&g, 10, 3, 2, &type_dist_mat, &pref_mat, &out_types, &in_types, IGRAPH_LOOPS);
    {
        /* Check that all out-types are 2 and all in-types are 0 */
        igraph_vector_int_t v;
        igraph_vector_int_init(&v, igraph_vcount(&g));

        igraph_vector_int_fill(&v, 2);
        IGRAPH_ASSERT(igraph_vector_int_all_e(&out_types, &v));

        igraph_vector_int_fill(&v, 0);
        IGRAPH_ASSERT(igraph_vector_int_all_e(&in_types, &v));

        igraph_vector_int_destroy(&v);
    }

    igraph_destroy(&g);

    igraph_matrix_destroy(&type_dist_mat);
    igraph_vector_destroy(&type_dist);
    igraph_matrix_destroy(&pref_mat);

    igraph_vector_int_destroy(&out_types);
    igraph_vector_int_destroy(&in_types);

    VERIFY_FINALLY_STACK();

    return 0;
}
