/* -*- mode: C -*-  */
/*
   IGraph library.
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

igraph_bool_t vector_eq(const igraph_vector_t *a, const igraph_vector_t *b) {
    igraph_integer_t na = igraph_vector_size(a);
    igraph_integer_t nb = igraph_vector_size(b);
    if (na != nb) {
        return false;
    }
    for (igraph_integer_t i=0; i < na; i++) {
        if (isnan(VECTOR(*a)[i]) && isnan(VECTOR(*b)[i])) {
            continue;
        }
        if (! igraph_almost_equals(VECTOR(*a)[i], VECTOR(*b)[i], 1e-12)) {
            return false;
        }
    }
    return true;
}

void check(const igraph_t *g, igraph_neimode_t mode1, igraph_neimode_t mode2, const igraph_vector_t *weights) {
    igraph_vector_t knn, knnk1, knnk2;

    igraph_vector_init(&knn, 0);
    igraph_vector_init(&knnk1, 0);
    igraph_vector_init(&knnk2, 0);

    igraph_avg_nearest_neighbor_degree(g, igraph_vss_all(),
                                       mode1, mode2,
                                       &knn, &knnk1, weights);

    printf("k_nn_i: ");
    print_vector(&knn);

    igraph_knnk(g, &knnk2, weights, mode1, mode2, mode1 == IGRAPH_ALL ? false : true);

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

int main(void) {
    igraph_t ug, dg;
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
    check(&ug, IGRAPH_ALL, IGRAPH_ALL, NULL);

    printf("Directed treated as undirected\n");
    check(&dg, IGRAPH_ALL, IGRAPH_ALL, NULL);

    printf("Directed: out, all\n");
    check(&dg, IGRAPH_OUT, IGRAPH_ALL, NULL);

    printf("Directed: out, in\n");
    check(&dg, IGRAPH_OUT, IGRAPH_IN, NULL);

    printf("Directed: out, out\n");
    check(&dg, IGRAPH_OUT, IGRAPH_OUT, NULL);

    printf("\nWEIGHTED\n\n");
    printf("Undirected\n");
    check(&ug, IGRAPH_ALL, IGRAPH_ALL, &uweights);

    printf("Directed treated as undirected\n");
    check(&dg, IGRAPH_ALL, IGRAPH_ALL, &dweights);

    printf("Directed: out, all\n");
    check(&dg, IGRAPH_OUT, IGRAPH_ALL, &dweights);

    printf("Directed: out, in\n");
    check(&dg, IGRAPH_OUT, IGRAPH_IN, &dweights);

    printf("Directed: out, out\n");
    check(&dg, IGRAPH_OUT, IGRAPH_OUT, &dweights);

    igraph_vector_destroy(&dweights);
    igraph_destroy(&dg);
    igraph_vector_destroy(&uweights);
    igraph_destroy(&ug);

    VERIFY_FINALLY_STACK();

    return 0;
}
