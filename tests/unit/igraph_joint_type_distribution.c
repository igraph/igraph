/*
   igraph library.
   Copyright (C) 2023  The igraph development team <igraph@igraph.org>

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

/* Compare with igraph_modularity() and igraph_assortativity_nominal() */
void check_assort(const igraph_t *g, const igraph_vector_t *weights, const igraph_vector_int_t *types) {
    igraph_vector_t a, b;
    igraph_matrix_t p;
    igraph_int_t n;
    igraph_real_t c1, c2, q1, q2;

    igraph_vector_init(&a, 0);
    igraph_vector_init(&b, 0);

    igraph_matrix_init(&p, 0, 0);

    igraph_joint_type_distribution(g, weights, &p, types, NULL, /*directed*/ true, /*normalized*/ true);

    igraph_matrix_rowsum(&p, &a);
    igraph_matrix_colsum(&p, &b);

    n = igraph_matrix_nrow(&p);
    IGRAPH_ASSERT(igraph_matrix_ncol(&p) == n);

    c1 = c2 = 0;
    for (igraph_int_t i=0; i < n; i++) {
        c1 += MATRIX(p, i, i);
        c2 += VECTOR(a)[i] * VECTOR(b)[i];
    }
    q2 = c1 - c2;

    igraph_modularity(g, types, weights, 1, true, &q1);
    // printf("Modularity: %g == %g\n", q1, q2);
    IGRAPH_ASSERT(igraph_almost_equals(q1, q2, 1e-14));

    if (! weights) {
        q1 /= 1 - c2;
        igraph_assortativity_nominal(g, NULL, types, &q2, /*directed*/ true, /*normalized*/ true);
        // printf("Normalized nominal assortativity: %g == %g\n", q1, q2);
        IGRAPH_ASSERT(igraph_almost_equals(q1, q2, 1e-14));
    }

    igraph_matrix_destroy(&p);

    igraph_vector_destroy(&b);
    igraph_vector_destroy(&a);
}

int main(void) {
    igraph_t graph;
    igraph_vector_int_t t1, t2;
    igraph_vector_t weights;
    igraph_matrix_t p;

    igraph_matrix_init(&p, 0, 0);

    printf("Null graph\n");
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    igraph_vector_int_init(&t1, 0);
    igraph_joint_type_distribution(&graph, NULL, &p, &t1, NULL, false, false);
    print_matrix(&p);
    igraph_vector_int_destroy(&t1);
    igraph_destroy(&graph);

    printf("\nSingleton graph\n");
    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    igraph_vector_int_init_int(&t1, 1, 0);
    igraph_joint_type_distribution(&graph, NULL, &p, &t1, NULL, false, false);
    print_matrix(&p);
    igraph_vector_int_destroy(&t1);
    igraph_destroy(&graph);

    printf("\nUndirected singleton with loop\n");
    igraph_small(&graph, 1, IGRAPH_UNDIRECTED, 0,0, -1);
    igraph_vector_int_init_int(&t1, 1, 0);
    igraph_joint_type_distribution(&graph, NULL, &p, &t1, NULL, false, false);
    print_matrix(&p);
    igraph_vector_int_destroy(&t1);
    igraph_destroy(&graph);

    printf("\nDirected singleton with loop\n");
    igraph_small(&graph, 1, IGRAPH_DIRECTED, 0,0, -1);
    igraph_vector_int_init_int(&t1, 1, 0);
    igraph_joint_type_distribution(&graph, NULL, &p, &t1, NULL, false, false);
    print_matrix(&p);
    igraph_vector_int_destroy(&t1);
    igraph_destroy(&graph);

    printf("\nSmall undirected graph\n");
    igraph_small(&graph, 6, IGRAPH_UNDIRECTED,
                 3, 0, 0, 3, 0, 2, 3, 1, 5, 5, 4, 2, 1, 1, 1, 1, 0, 1, 5, 1,
                 -1);

    igraph_vector_int_init_int(&t1, 6,
                               0, 0, 1, 1, 2, 2);

    igraph_joint_type_distribution(&graph, NULL, &p, &t1, NULL, false, false);
    print_matrix(&p);

    IGRAPH_ASSERT(igraph_matrix_sum(&p) == 2*igraph_ecount(&graph));
    check_assort(&graph, NULL, &t1);

    igraph_vector_init_range(&weights, 10, 10 + igraph_ecount(&graph));
    check_assort(&graph, &weights, &t1);
    igraph_vector_destroy(&weights);

    igraph_vector_int_destroy(&t1);
    igraph_destroy(&graph);

    printf("\nSmall directed graph\n");
    igraph_small(&graph, 6, IGRAPH_DIRECTED,
                 1, 1, 2, 4, 0, 2, 2, 2, 3, 2, 3, 0, 2, 1, 2, 3, 4, 1, 2, 4,
                 -1);

    igraph_vector_int_init_int(&t1, 6,
                               0, 0, 1, 1, 2, 2);

    igraph_joint_type_distribution(&graph, NULL, &p, &t1, NULL, true, false);

    printf("From and to types are the same:\n");
    print_matrix(&p);
    IGRAPH_ASSERT(igraph_matrix_sum(&p) == igraph_ecount(&graph));

    check_assort(&graph, NULL, &t1);

    igraph_vector_init_range(&weights, 10, 10 + igraph_ecount(&graph));
    check_assort(&graph, &weights, &t1);
    igraph_vector_destroy(&weights);

    igraph_vector_int_init_int(&t2, 6,
                               0, 1, 1, 1, 0, 1);
    printf("From and to types are different:\n");
    igraph_joint_type_distribution(&graph, NULL, &p, &t1, &t2, true, false);
    print_matrix(&p);
    IGRAPH_ASSERT(igraph_matrix_sum(&p) == igraph_ecount(&graph));

    igraph_vector_int_destroy(&t2);
    igraph_vector_int_destroy(&t1);
    igraph_destroy(&graph);

    igraph_matrix_destroy(&p);

    VERIFY_FINALLY_STACK();
    return 0;
}
