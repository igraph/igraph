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

igraph_bool_t check_laplacian(igraph_t* graph, const igraph_matrix_t* matrix, const igraph_vector_t* w) {
    igraph_vector_t vec, res;
    long int i, j;

    igraph_vector_init(&vec, 0);
    igraph_vector_init(&res, igraph_vcount(graph));

    if (w) {
        igraph_strength(graph, &vec, igraph_vss_all(), IGRAPH_OUT, IGRAPH_NO_LOOPS, w);
    } else {
        igraph_degree(graph, &vec, igraph_vss_all(), IGRAPH_OUT, IGRAPH_NO_LOOPS);
    }

    for (i = 0; i < igraph_vcount(graph); i++) {
        VECTOR(vec)[i] = sqrt(VECTOR(vec)[i]);
    }

    for (i = 0; i < igraph_vcount(graph); i++) {
        for (j = 0; j < igraph_vcount(graph); j++) {
            VECTOR(res)[i] += MATRIX(*matrix, i, j) * VECTOR(vec)[j];
        }
    }

    if (igraph_vector_min(&res) > 1e-7) {
        printf("Invalid Laplacian matrix:\n");
        igraph_matrix_print(matrix);
        return 0;
    }

    igraph_vector_destroy(&vec);
    igraph_vector_destroy(&res);

    return 1;
}

int test_unnormalized_laplacian(const igraph_vector_t* w, igraph_bool_t dir) {
    igraph_t g;
    igraph_matrix_t m, m2;
    igraph_sparsemat_t sm;
    igraph_vector_t vec, *weights = NULL;
    igraph_matrix_init(&m, 1, 1);
    igraph_sparsemat_init(&sm, 0, 0, 0);

    if (w) {
        weights = (igraph_vector_t*)calloc(1, sizeof(igraph_vector_t));
        igraph_vector_copy(weights, w);
    }

    /* No loop or multiple edges */
    igraph_ring(&g, 5, dir, 0, 1);
    igraph_laplacian(&g, &m, &sm, 0, weights);
    igraph_matrix_init(&m2, 0, 0);
    igraph_sparsemat_as_matrix(&m2, &sm);
    if (!igraph_matrix_all_e_tol(&m, &m2, 0)) {
        return 41;
    }
    igraph_matrix_destroy(&m2);
    igraph_matrix_print(&m);
    printf("===\n");

    /* Add some loop edges */
    igraph_vector_init_real(&vec, 4, 1.0, 1.0, 2.0, 2.0);
    igraph_add_edges(&g, &vec, 0);
    igraph_vector_destroy(&vec);
    if (weights) {
        igraph_vector_push_back(weights, 2);
        igraph_vector_push_back(weights, 2);
    }

    igraph_laplacian(&g, &m, &sm, 0, weights);
    igraph_matrix_init(&m2, 0, 0);
    igraph_sparsemat_as_matrix(&m2, &sm);
    if (!igraph_matrix_all_e_tol(&m, &m2, 0)) {
        return 42;
    }
    igraph_matrix_destroy(&m2);
    igraph_matrix_print(&m);
    printf("===\n");

    /* Duplicate some edges */
    igraph_vector_init_real(&vec, 4, 1.0, 2.0, 3.0, 4.0);
    igraph_add_edges(&g, &vec, 0);
    igraph_vector_destroy(&vec);
    if (weights) {
        igraph_vector_push_back(weights, 3);
        igraph_vector_push_back(weights, 3);
    }

    igraph_laplacian(&g, &m, &sm, 0, weights);
    igraph_matrix_init(&m2, 0, 0);
    igraph_sparsemat_as_matrix(&m2, &sm);
    if (!igraph_matrix_all_e_tol(&m, &m2, 0)) {
        return 43;
    }
    igraph_matrix_destroy(&m2);
    igraph_matrix_print(&m);

    igraph_destroy(&g);

    igraph_matrix_destroy(&m);
    if (weights) {
        igraph_vector_destroy(weights);
        free(weights);
    }

    igraph_sparsemat_destroy(&sm);

    return 0;
}

int test_normalized_laplacian(const igraph_vector_t *w, igraph_bool_t dir) {
    igraph_t g;
    igraph_matrix_t m, m2;
    igraph_sparsemat_t sm;
    igraph_vector_t vec, *weights = 0;
    igraph_bool_t ok = 1;
    igraph_matrix_init(&m, 1, 1);
    igraph_sparsemat_init(&sm, 0, 0, 0);

    if (w) {
        weights = (igraph_vector_t*) calloc(1, sizeof(igraph_vector_t));
        igraph_vector_copy(weights, w);
    }

    /* Undirected graph, no loop or multiple edges */
    igraph_ring(&g, 5, dir, 0, 1);
    igraph_laplacian(&g, &m, &sm, 1, weights);
    igraph_matrix_init(&m2, 0, 0);
    igraph_sparsemat_as_matrix(&m2, &sm);
    if (!igraph_matrix_all_e_tol(&m, &m2, 0)) {
        return 44;
    }
    igraph_matrix_destroy(&m2);
    ok = ok && check_laplacian(&g, &m, weights);

    /* Add some loop edges */
    igraph_vector_init_real(&vec, 4, 1.0, 1.0, 2.0, 2.0);
    igraph_add_edges(&g, &vec, 0);
    igraph_vector_destroy(&vec);
    if (weights) {
        igraph_vector_push_back(weights, 2);
        igraph_vector_push_back(weights, 2);
    }

    igraph_laplacian(&g, &m, &sm, 1, weights);
    igraph_matrix_init(&m2, 0, 0);
    igraph_sparsemat_as_matrix(&m2, &sm);
    if (!igraph_matrix_all_e_tol(&m, &m2, 0)) {
        return 45;
    }
    igraph_matrix_destroy(&m2);
    ok = ok && check_laplacian(&g, &m, weights);

    /* Duplicate some edges */
    igraph_vector_init_real(&vec, 4, 1.0, 2.0, 3.0, 4.0);
    igraph_add_edges(&g, &vec, 0);
    igraph_vector_destroy(&vec);
    if (weights) {
        igraph_vector_push_back(weights, 3);
        igraph_vector_push_back(weights, 3);
    }

    igraph_laplacian(&g, &m, &sm, 1, weights);
    igraph_matrix_init(&m2, 0, 0);
    igraph_sparsemat_as_matrix(&m2, &sm);
    if (!igraph_matrix_all_e_tol(&m, &m2, 0)) {
        return 46;
    }
    igraph_matrix_destroy(&m2);
    ok = ok && check_laplacian(&g, &m, weights);

    igraph_destroy(&g);

    igraph_matrix_destroy(&m);
    if (weights) {
        igraph_vector_destroy(weights);
        free(weights);
    }

    if (ok) {
        printf("OK\n");
    }

    igraph_sparsemat_destroy(&sm);

    return !ok;
}

int main() {
    int res;
    int i;
    igraph_vector_t weights;

    igraph_vector_init_real(&weights, 5, 1.0, 2.0, 3.0, 4.0, 5.0);

    for (i = 0; i < 8; i++) {
        igraph_bool_t is_normalized = i / 4;
        igraph_vector_t* v = ((i & 2) / 2 ? &weights : 0);
        igraph_bool_t dir = (i % 2 ? IGRAPH_DIRECTED : IGRAPH_UNDIRECTED);

        printf("=== %sormalized, %sweighted, %sdirected\n",
               (is_normalized ? "N" : "Unn"),
               (v != 0 ? "" : "un"),
               (dir == IGRAPH_DIRECTED ? "" : "un")
              );

        if (is_normalized) {
            res = test_normalized_laplacian(v, dir);
        } else {
            res = test_unnormalized_laplacian(v, dir);
        }

        if (res) {
            return i + 1;
        }
    }

    igraph_vector_destroy(&weights);

    return 0;
}
