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

#include <igraph.h>
#include "test_utilities.h"

/* Compare vector for equality, allowing for NaN elements.
 * Note that NaN != NaN, thus we cannot use vector_all_e(). */
igraph_bool_t vec_equal(const igraph_vector_t *v1, const igraph_vector_t *v2) {
    igraph_integer_t n = igraph_vector_size(v1);

    if (igraph_vector_size(v2) != n) return false;

    for (igraph_integer_t i=0; i < n; i++) {
        igraph_real_t x1 = VECTOR(*v1)[i], x2 = VECTOR(*v2)[i];

        if (isnan(x1)) {
            if (! isnan(x2)) {
                return false;
            }
        } else if (x1 != x2) {
            return false;
        }
    }
    return true;
}

/* This is a trivial implementation of ECC for k=3 using igraph_list_triangles() */
void get_ecc3(const igraph_t *g, igraph_vector_t *res, igraph_bool_t offset, igraph_bool_t normalize) {
    igraph_vector_int_t triangles_vec, eids;
    igraph_matrix_int_t triangles;

    igraph_vector_resize(res, igraph_ecount(g));
    igraph_vector_null(res);

    igraph_vector_int_init(&triangles_vec, 0);
    igraph_list_triangles(g, &triangles_vec);

    igraph_matrix_int_view_from_vector(&triangles, &triangles_vec, 3);

    igraph_vector_int_init(&eids, 0);

    for (igraph_integer_t i=0; i < igraph_matrix_int_ncol(&triangles); i++) {
        igraph_integer_t u, v, w;
        igraph_integer_t ec;

        u = MATRIX(triangles, 0, i);
        v = MATRIX(triangles, 1, i);
        w = MATRIX(triangles, 2, i);

        igraph_get_all_eids_between(g, &eids, u, v, IGRAPH_UNDIRECTED);
        ec = igraph_vector_int_size(&eids);
        for (igraph_integer_t j=0; j < ec; j++) {
            VECTOR(*res)[ VECTOR(eids)[j] ] += 1;
        }

        igraph_get_all_eids_between(g, &eids, v, w, IGRAPH_UNDIRECTED);
        ec = igraph_vector_int_size(&eids);
        for (igraph_integer_t j=0; j < ec; j++) {
            VECTOR(*res)[ VECTOR(eids)[j] ] += 1;
        }

        igraph_get_all_eids_between(g, &eids, w, u, IGRAPH_UNDIRECTED);
        ec = igraph_vector_int_size(&eids);
        for (igraph_integer_t j=0; j < ec; j++) {
            VECTOR(*res)[ VECTOR(eids)[j] ] += 1;
        }
    }

    igraph_vector_int_t degree;
    igraph_vector_int_init(&degree, 0);
    igraph_degree(g, &degree, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);

    for (igraph_integer_t e=0; e < igraph_ecount(g); e++) {
        igraph_integer_t v1 = IGRAPH_FROM(g, e), v2 = IGRAPH_TO(g, e);
        igraph_real_t s;
        if (v1 == v2) {
            s = 0.0;
        } else {
            igraph_integer_t d1 = VECTOR(degree)[v1], d2 = VECTOR(degree)[v2];
            s = (d1 < d2 ? d1 : d2) - 1.0;
        }
        if (offset) VECTOR(*res)[e] += 1;
        if (normalize) VECTOR(*res)[e] /= s;
    }

    igraph_vector_int_destroy(&degree);

    igraph_vector_int_destroy(&eids);
    igraph_vector_int_destroy(&triangles_vec);
}

/* These are globals, as they're used both from test_ecc() and main() */
igraph_vector_t ecc1, ecc2, ecc3;

void test_ecc(const igraph_t *g) {
    printf("k=3\n");
    igraph_ecc(g, &ecc1, igraph_ess_all(IGRAPH_EDGEORDER_ID), 3, false, true);
    print_vector(&ecc1);

    igraph_ecc(g, &ecc2, igraph_ess_range(0, igraph_ecount(g)), 3, false, true);
    print_vector(&ecc2);
    IGRAPH_ASSERT(vec_equal(&ecc1, &ecc2));

    get_ecc3(g, &ecc3, false, true);
    print_vector(&ecc3);

    IGRAPH_ASSERT(vec_equal(&ecc1, &ecc3));

    printf("\nk=4\n");
    igraph_ecc(g, &ecc1, igraph_ess_all(IGRAPH_EDGEORDER_ID), 4, false, true);
    print_vector(&ecc1);

    igraph_ecc(g, &ecc2, igraph_ess_range(0, igraph_ecount(g)), 4, false, true);
    print_vector(&ecc2);
    IGRAPH_ASSERT(vec_equal(&ecc1, &ecc2));
}

int main(void) {
    igraph_t g;

    igraph_vector_init(&ecc1, 0);
    igraph_vector_init(&ecc2, 0);
    igraph_vector_init(&ecc3, 0);

    /* Null graph */

    printf("Null graph:\n");
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    test_ecc(&g);
    igraph_destroy(&g);

    /* Singleton */

    printf("\nSingleton graph:\n");
    igraph_empty(&g, 1, IGRAPH_UNDIRECTED);
    test_ecc(&g);
    igraph_destroy(&g);

    /* P_2 */

    printf("\nP_2 graph\n");
    igraph_full(&g, 2, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    test_ecc(&g);
    igraph_destroy(&g);

    /* K_5 */

    printf("\nK_5 graph:\n");
    igraph_full(&g, 5, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    test_ecc(&g);
    igraph_destroy(&g);

    printf("\nK_5 graph with self-loops:\n");
    igraph_full(&g, 5, IGRAPH_UNDIRECTED, IGRAPH_LOOPS);
    test_ecc(&g);
    igraph_destroy(&g);

    /* Multigraph */

    printf("\nMultigraph:\n");
    igraph_small(&g, 5, IGRAPH_UNDIRECTED,
        0,1, 1,2, 2,0, 0,1, 1,3, 3,4, 4,0, 0,5, 5,5, 5,5, 1,4,
        -1);
    test_ecc(&g);
    igraph_destroy(&g);

    /* Karate club */

    printf("\nZachary karate club:\n");
    igraph_famous(&g, "Zachary");

    test_ecc(&g);

    /* Check invalid input */

    CHECK_ERROR(igraph_ecc(&g, &ecc1, igraph_ess_all(IGRAPH_EDGEORDER_ID), 2, false, true), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_ecc(&g, &ecc1, igraph_ess_all(IGRAPH_EDGEORDER_ID), 5, false, true), IGRAPH_UNIMPLEMENTED);

    igraph_destroy(&g);

    igraph_vector_destroy(&ecc3);
    igraph_vector_destroy(&ecc2);
    igraph_vector_destroy(&ecc1);

    VERIFY_FINALLY_STACK();

    return 0;
}
