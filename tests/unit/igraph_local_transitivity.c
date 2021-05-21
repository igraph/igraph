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
#include <math.h>

#include "test_utilities.inc"

/* Compare the elements of two vectors for equality, handling NaN values. */
igraph_bool_t vector_equal(const igraph_vector_t *v1, const igraph_vector_t *v2) {
    long int n1 = igraph_vector_size(v1), n2 = igraph_vector_size(v2);
    long int i;

    if (n1 != n2) {
        return 0;
    }

    for (i=0; i < n1; ++i) {
        /* Since NaN == NaN compares false, we must handle NaN values early. */
        if (igraph_is_nan(VECTOR(*v1)[i]) && igraph_is_nan(VECTOR(*v2)[i])) {
            continue;
        }
        if (VECTOR(*v1)[i]  != VECTOR(*v2)[i]) {
            return 0;
        }
    }

    return 1;
}

/* Compute the average of a vector, ignoring NaN values. */
igraph_real_t vector_avg(const igraph_vector_t *v) {
    long int n = igraph_vector_size(v);
    long int i;
    igraph_real_t sum = 0.0, count;

    count = 0;
    for (i=0; i < n; ++i) {
        if (igraph_is_nan(VECTOR(*v)[i])) {
            continue;
        }
        sum += VECTOR(*v)[i];
        count += 1;
    }
    return sum / count;
}

int main() {

    igraph_t g;
    igraph_vector_t result1, result2, result3;
    igraph_vs_t vertices;
    igraph_real_t avg_local;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_vector_init(&result1, 0);
    igraph_vector_init(&result2, 0);

    /* igraph_transitivity_local_undirected() uses different code paths for:
     *  - all vertices
     *  - some vertices of graphs with >= 100 vertices
     *  - some vertices of graphs with < 100 vertices
     *
     * We test that these are consistent.
     */

    /* 100 vertices */

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 100, 0.1,
                            IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

    igraph_vs_seq(&vertices, 0, igraph_vcount(&g) - 1);

    igraph_transitivity_local_undirected(&g, &result1, igraph_vss_all(),
                                         IGRAPH_TRANSITIVITY_NAN);
    igraph_transitivity_local_undirected(&g, &result2, vertices,
                                         IGRAPH_TRANSITIVITY_NAN);

    IGRAPH_ASSERT(vector_equal(&result1, &result2));

    igraph_vs_destroy(&vertices);
    igraph_destroy(&g);

    /* 50 vertices */

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 50, 0.3,
                            IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

    igraph_vs_seq(&vertices, 0, igraph_vcount(&g) - 1);

    igraph_transitivity_local_undirected(&g, &result1, igraph_vss_all(),
                                         IGRAPH_TRANSITIVITY_NAN);
    igraph_transitivity_local_undirected(&g, &result2, vertices,
                                         IGRAPH_TRANSITIVITY_NAN);

    IGRAPH_ASSERT(vector_equal(&result1, &result2));

    igraph_vs_destroy(&vertices);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    /* Zachary karate club */

    printf("Zachary karate club network:\n");
    igraph_famous(&g, "Zachary");

    printf("IGRAPH_TRANSITIVITY_ZERO:\n");
    igraph_transitivity_local_undirected(&g, &result1, igraph_vss_all(), IGRAPH_TRANSITIVITY_ZERO);
    print_vector(&result1);

    printf("IGRAPH_TRANSITIVITY_NAN:\n");
    igraph_transitivity_local_undirected(&g, &result1, igraph_vss_all(), IGRAPH_TRANSITIVITY_NAN);
    print_vector(&result1);

    igraph_destroy(&g);

    /* Small graphs */

    printf("\nNull graph:\n");
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    igraph_transitivity_local_undirected(&g, &result1, igraph_vss_all(), IGRAPH_TRANSITIVITY_NAN);
    print_vector(&result1);
    igraph_transitivity_avglocal_undirected(&g, &avg_local, IGRAPH_TRANSITIVITY_NAN);
    IGRAPH_ASSERT(igraph_is_nan(avg_local));
    igraph_transitivity_avglocal_undirected(&g, &avg_local, IGRAPH_TRANSITIVITY_ZERO);
    IGRAPH_ASSERT(avg_local == 0);
    igraph_destroy(&g);

    printf("\nSingleton graph:\n");
    igraph_empty(&g, 1, IGRAPH_UNDIRECTED);
    igraph_transitivity_local_undirected(&g, &result1, igraph_vss_all(), IGRAPH_TRANSITIVITY_NAN);
    print_vector(&result1);
    igraph_transitivity_avglocal_undirected(&g, &avg_local, IGRAPH_TRANSITIVITY_NAN);
    IGRAPH_ASSERT(igraph_is_nan(avg_local));
    igraph_transitivity_avglocal_undirected(&g, &avg_local, IGRAPH_TRANSITIVITY_ZERO);
    IGRAPH_ASSERT(avg_local == 0);
    igraph_destroy(&g);

    printf("\nTwo connected vertices:\n");
    igraph_small(&g, 2, IGRAPH_UNDIRECTED, 0,1, -1);
    igraph_transitivity_local_undirected(&g, &result1, igraph_vss_all(), IGRAPH_TRANSITIVITY_NAN);
    print_vector(&result1);
    igraph_destroy(&g);

    printf("\nTriangle:\n");
    igraph_full(&g, 3, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_transitivity_local_undirected(&g, &result1, igraph_vss_all(), IGRAPH_TRANSITIVITY_NAN);
    print_vector(&result1);
    igraph_destroy(&g);

    printf("\nTwo-star:\n");
    igraph_small(&g, 3, IGRAPH_UNDIRECTED, 0,2, 0,1, -1);
    igraph_transitivity_local_undirected(&g, &result1, igraph_vss_all(), IGRAPH_TRANSITIVITY_NAN);
    print_vector(&result1);
    igraph_destroy(&g);

    /* Multigraph */

    printf("\nDirected and multigraphs:\n");

    igraph_small(&g, 20, IGRAPH_DIRECTED,
                 15, 12, 12, 10, 15, 0, 11, 10, 2, 8, 8, 6, 13, 17, 10, 10, 17, 2, 14,
                 0, 16, 13, 14, 14, 0, 5, 6, 4, 0, 9, 0, 6, 10, 9, 16, 4, 14, 5, 17,
                 15, 14, 9, 17, 17, 1, 4, 10, 16, 7, 0, 11, 12, 6, 13, 2, 17, 4, 0, 0,
                 14, 4, 0, 6, 16, 16, 14, 13, 13, 12, 11, 3, 11, 11, 3, 6, 7, 4, 14,
                 10, 8, 13, 7, 14, 2, 5, 2, 0, 14, 3, 15, 5, 5, 7, 2, 14, 15, 5, 10,
                 10, 16, 7, 9, 14, 0, 15, 7, 13, 1, 15, 1, 4, 5, 4, 6, 16, 13, 6, 17,
                 8, 6, 9, 3, 8, 6, 6, 14, 11, 14, 6, 10, 10, 5, 1, 0, 16, 17, 9, 1, 5,
                 0, 5, 15, 8, 0, 0, 8, 5, 3, 9, 4, 13, 12, 11, 0, 11, 0, 10, 6, 4, 13,
                 8, 9, 11, 11, 3, 16, 1, 2, 16, 0, 9, 8, 3, 8, 8, 7, 12, 10, 9, 3, 13,
                 5, 3, 9, 6, 2, 11, 10, 1, 16, 0, 2, 10, 17, 16, 8, 11, 5, 13, 0, 19, 19,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 -1);

    igraph_vs_seq(&vertices, 0, igraph_vcount(&g) - 1);

    printf("\nDirected multi:\n");
    igraph_transitivity_local_undirected(&g, &result1, igraph_vss_all(), IGRAPH_TRANSITIVITY_NAN);
    print_vector(&result1);

    igraph_vector_copy(&result3, &result1);

    igraph_transitivity_local_undirected(&g, &result2, vertices, IGRAPH_TRANSITIVITY_NAN);
    print_vector(&result2);
    IGRAPH_ASSERT(vector_equal(&result2, &result3));

    igraph_transitivity_avglocal_undirected(&g, &avg_local, IGRAPH_TRANSITIVITY_NAN);
    printf("Average: %.10g == %.10g == %.10g\n", avg_local, vector_avg(&result1), vector_avg(&result2));
    IGRAPH_ASSERT(fabs(avg_local - vector_avg(&result1)) < 1e-14);

    printf("\nUndirected multi:\n");
    igraph_to_undirected(&g, IGRAPH_TO_UNDIRECTED_COLLAPSE, NULL);
    igraph_transitivity_local_undirected(&g, &result1, igraph_vss_all(), IGRAPH_TRANSITIVITY_NAN);
    print_vector(&result1);
    IGRAPH_ASSERT(vector_equal(&result1, &result3));
    igraph_transitivity_local_undirected(&g, &result2, vertices, IGRAPH_TRANSITIVITY_NAN);
    print_vector(&result2);
    IGRAPH_ASSERT(vector_equal(&result2, &result3));

    igraph_transitivity_avglocal_undirected(&g, &avg_local, IGRAPH_TRANSITIVITY_NAN);
    printf("Average: %.10g == %.10g == %.10g\n", avg_local, vector_avg(&result1), vector_avg(&result2));
    IGRAPH_ASSERT(fabs(avg_local - vector_avg(&result1)) < 1e-14);

    printf("\nUndirected simple:\n");
    igraph_simplify(&g, 1, 1, NULL);
    igraph_transitivity_local_undirected(&g, &result1, igraph_vss_all(), IGRAPH_TRANSITIVITY_NAN);
    print_vector(&result1);
    IGRAPH_ASSERT(vector_equal(&result1, &result3));
    igraph_transitivity_local_undirected(&g, &result2, vertices, IGRAPH_TRANSITIVITY_NAN);
    print_vector(&result2);
    IGRAPH_ASSERT(vector_equal(&result2, &result3));

    igraph_transitivity_avglocal_undirected(&g, &avg_local, IGRAPH_TRANSITIVITY_NAN);
    printf("Average: %.10g == %.10g == %.10g\n", avg_local, vector_avg(&result1), vector_avg(&result2));
    IGRAPH_ASSERT(fabs(avg_local - vector_avg(&result1)) < 1e-14);

    igraph_vector_destroy(&result3);
    igraph_vector_destroy(&result2);
    igraph_vector_destroy(&result1);

    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
