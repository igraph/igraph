/*
   igraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

/* Check vector equality with tolerances. Consider NaN values equal. */
igraph_bool_t vector_eq(const igraph_vector_t *a, const igraph_vector_t *b) {
    igraph_int_t na = igraph_vector_size(a);
    igraph_int_t nb = igraph_vector_size(b);
    if (na != nb) {
        return false;
    }
    for (igraph_int_t i=0; i < na; i++) {
        if (isnan(VECTOR(*a)[i]) && isnan(VECTOR(*b)[i])) {
            continue;
        }
        if (! igraph_almost_equals(VECTOR(*a)[i], VECTOR(*b)[i], 1e-12)) {
            return false;
        }
    }
    return true;
}

/* Compare results from igraph_avg_nearest_neighbor_degree() and igraph_degree_correlation_vector() */
void compare_implementations(const igraph_t *g, igraph_neimode_t mode1, igraph_neimode_t mode2, const igraph_vector_t *weights) {
    igraph_vector_t knn, knnk1, knnk2;

    igraph_vector_init(&knn, 0);
    igraph_vector_init(&knnk1, 0);
    igraph_vector_init(&knnk2, 0);

    igraph_avg_nearest_neighbor_degree(g, igraph_vss_all(),
                                       mode1, mode2,
                                       &knn, &knnk1, weights);

    printf("k_nn_i: ");
    print_vector(&knn);

    igraph_degree_correlation_vector(g, weights, &knnk2, mode1, mode2, mode1 == IGRAPH_ALL ? false : true);

    printf("k_nn(k): ");
    print_vector(&knnk2);
    /* print_vector(&knnk1); */

    IGRAPH_ASSERT(isnan(VECTOR(knnk2)[0]));

    igraph_vector_remove(&knnk2, 0);
    IGRAPH_ASSERT(vector_eq(&knnk1, &knnk2));

    igraph_vector_destroy(&knnk2);
    igraph_vector_destroy(&knnk1);
    igraph_vector_destroy(&knn);

    printf("\n");
}

/* Check that undirected graphs are treated as a directed ones with all-reciprocal edges.
 * This helps verify non-simple graphs, especially those with self-loops. */
void check_dir_undir_equiv(const igraph_t *ug) {
    igraph_t rg;
    igraph_vector_t knn1, knn2, knnk1, knnk2;

    igraph_vector_init(&knn1, 0);
    igraph_vector_init(&knn2, 0);
    igraph_vector_init(&knnk1, 0);
    igraph_vector_init(&knnk2, 0);

    igraph_copy(&rg, ug);
    igraph_to_directed(&rg, IGRAPH_TO_DIRECTED_MUTUAL);

    igraph_avg_nearest_neighbor_degree(ug, igraph_vss_all(), IGRAPH_ALL, IGRAPH_ALL, &knn1, &knnk1, NULL);
    igraph_avg_nearest_neighbor_degree(&rg, igraph_vss_all(), IGRAPH_OUT, IGRAPH_OUT, &knn2, &knnk2, NULL);

    IGRAPH_ASSERT(vector_eq(&knn1, &knn2));
    IGRAPH_ASSERT(vector_eq(&knnk1, &knnk2));

    igraph_avg_nearest_neighbor_degree(&rg, igraph_vss_all(), IGRAPH_IN, IGRAPH_IN, &knn2, &knnk2, NULL);

    IGRAPH_ASSERT(vector_eq(&knn1, &knn2));
    IGRAPH_ASSERT(vector_eq(&knnk1, &knnk2));

    igraph_vector_null(&knnk1);
    igraph_vector_null(&knnk2);

    igraph_degree_correlation_vector(ug, NULL, &knnk1, IGRAPH_ALL, IGRAPH_ALL, false);
    igraph_degree_correlation_vector(&rg, NULL, &knnk2, IGRAPH_OUT, IGRAPH_OUT, false);

    IGRAPH_ASSERT(vector_eq(&knnk1, &knnk2));

    igraph_destroy(&rg);

    igraph_vector_destroy(&knnk2);
    igraph_vector_destroy(&knnk1);
    igraph_vector_destroy(&knn2);
    igraph_vector_destroy(&knn1);
}

int main(void) {
    igraph_t ug, dg, g;
    igraph_vector_t uweights, dweights;

    igraph_small(&ug, 10, IGRAPH_UNDIRECTED,
                 0, 4, 1, 1, 1, 2, 1, 2, 1, 8, 2, 3, 2, 6, 3,
                 5, 3, 6, 3, 7, 4, 7, 7, 8, 7, 8, 7, 8, 6, 6,
                 -1);

    igraph_vector_init_range(&uweights, 0, igraph_ecount(&ug));
    igraph_vector_scale(&uweights, 1.0/8);

    igraph_small(&dg, 6, IGRAPH_DIRECTED,
                 0, 0, 0, 1, 0, 2, 1, 0, 1, 3, 1, 4, 2, 0, 2,
                 3, 3, 2, 3, 3, 4, 2, 4, 3, 0, 1, 1, 3, 1, 3,
                 -1);

    igraph_vector_init_range(&dweights, 0, igraph_ecount(&dg));
    igraph_vector_scale(&dweights, 1.0/8);

    printf("UNWEIGHTED\n\n");
    printf("Undirected\n");
    compare_implementations(&ug, IGRAPH_ALL, IGRAPH_ALL, NULL);

    printf("Directed treated as undirected\n");
    compare_implementations(&dg, IGRAPH_ALL, IGRAPH_ALL, NULL);

    printf("Directed: out, all\n");
    compare_implementations(&dg, IGRAPH_OUT, IGRAPH_ALL, NULL);

    printf("Directed: out, in\n");
    compare_implementations(&dg, IGRAPH_OUT, IGRAPH_IN, NULL);

    printf("Directed: out, out\n");
    compare_implementations(&dg, IGRAPH_OUT, IGRAPH_OUT, NULL);

    printf("\nWEIGHTED\n\n");
    printf("Undirected\n");
    compare_implementations(&ug, IGRAPH_ALL, IGRAPH_ALL, &uweights);

    printf("Directed treated as undirected\n");
    compare_implementations(&dg, IGRAPH_ALL, IGRAPH_ALL, &dweights);

    printf("Directed: out, all\n");
    compare_implementations(&dg, IGRAPH_OUT, IGRAPH_ALL, &dweights);

    printf("Directed: out, in\n");
    compare_implementations(&dg, IGRAPH_OUT, IGRAPH_IN, &dweights);

    printf("Directed: out, out\n");
    compare_implementations(&dg, IGRAPH_OUT, IGRAPH_OUT, &dweights);

    printf("\nSPECIAL CASES\n\n");
    printf("Null graph\n");
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    compare_implementations(&g, IGRAPH_ALL, IGRAPH_ALL, NULL);
    igraph_destroy(&g);

    printf("Singleton with loop\n");
    igraph_small(&g, 1, IGRAPH_UNDIRECTED, 0,0, -1);
    compare_implementations(&g, IGRAPH_ALL, IGRAPH_ALL, NULL);
    igraph_destroy(&g);

    printf("Two isolated vertices\n");
    igraph_empty(&g, 2, IGRAPH_UNDIRECTED);
    compare_implementations(&g, IGRAPH_ALL, IGRAPH_ALL, NULL);
    igraph_destroy(&g);

    printf("Check directed/undirected equivalence principle\n");
    check_dir_undir_equiv(&ug);

    igraph_vector_destroy(&dweights);
    igraph_destroy(&dg);
    igraph_vector_destroy(&uweights);
    igraph_destroy(&ug);

    VERIFY_FINALLY_STACK();

    return 0;
}
