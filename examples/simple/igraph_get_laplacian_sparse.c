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

int test_laplacian(const igraph_vector_t *w, igraph_bool_t dir, igraph_laplacian_normalization_t normalization) {
    igraph_t g;
    igraph_matrix_t m;
    igraph_sparsemat_t sm;
    igraph_vector_int_t vec;
    igraph_vector_t *weights = 0;
    igraph_neimode_t mode = IGRAPH_OUT;

    igraph_sparsemat_init(&sm, 0, 0, 0);

    if (w) {
        weights = (igraph_vector_t*) calloc(1, sizeof(igraph_vector_t));
        igraph_vector_init_copy(weights, w);
    }

    /* Base graph, no loop or multiple edges */
    igraph_ring(&g, 5, dir, 0, 1);
    igraph_get_laplacian_sparse(&g, &sm, mode, normalization, weights);
    igraph_matrix_init(&m, 0, 0);
    igraph_sparsemat_as_matrix(&m, &sm);
    igraph_matrix_print(&m);
    igraph_matrix_destroy(&m);
    printf("===\n");

    /* Add some loop edges */
    igraph_vector_int_init_int(&vec, 4, 1, 1, 2, 2);
    igraph_add_edges(&g, &vec, 0);
    igraph_vector_int_destroy(&vec);
    if (weights) {
        igraph_vector_push_back(weights, 2);
        igraph_vector_push_back(weights, 2);
    }

    igraph_get_laplacian_sparse(&g, &sm, mode, normalization, weights);
    igraph_matrix_init(&m, 0, 0);
    igraph_sparsemat_as_matrix(&m, &sm);
    igraph_matrix_print(&m);
    igraph_matrix_destroy(&m);
    printf("===\n");

    /* Duplicate some edges */
    igraph_vector_int_init_int(&vec, 4, 1, 2, 3, 4);
    igraph_add_edges(&g, &vec, 0);
    igraph_vector_int_destroy(&vec);
    if (weights) {
        igraph_vector_push_back(weights, 3);
        igraph_vector_push_back(weights, 3);
    }

    igraph_get_laplacian_sparse(&g, &sm, mode, normalization, weights);
    igraph_matrix_init(&m, 0, 0);
    igraph_sparsemat_as_matrix(&m, &sm);
    igraph_matrix_print(&m);
    igraph_matrix_destroy(&m);
    printf("===\n");

    /* Add an isolated vertex */
    igraph_add_vertices(&g, 1, NULL);

    igraph_get_laplacian_sparse(&g, &sm, mode, normalization, weights);
    igraph_matrix_init(&m, 0, 0);
    igraph_sparsemat_as_matrix(&m, &sm);
    igraph_matrix_print(&m);
    igraph_matrix_destroy(&m);

    igraph_destroy(&g);

    if (weights) {
        igraph_vector_destroy(weights);
        free(weights);
    }

    igraph_sparsemat_destroy(&sm);

    return 0;
}

int main(void) {
    int res;
    int i;
    igraph_vector_t weights;

    igraph_vector_init_int(&weights, 5, 1, 2, 3, 4, 5);

    for (i = 0; i < 8; i++) {
        igraph_bool_t is_normalized = i / 4;
        igraph_vector_t* v = ((i & 2) / 2 ? &weights : 0);
        igraph_bool_t dir = (i % 2 ? IGRAPH_DIRECTED : IGRAPH_UNDIRECTED);

        printf("=== %sormalized, %sweighted, %sdirected\n",
               (is_normalized ? "N" : "Unn"),
               (v != 0 ? "" : "un"),
               (dir == IGRAPH_DIRECTED ? "" : "un")
              );

        res = test_laplacian(v, dir, is_normalized ? IGRAPH_LAPLACIAN_SYMMETRIC : IGRAPH_LAPLACIAN_UNNORMALIZED);

        if (res) {
            return i + 1;
        }
    }

    igraph_vector_destroy(&weights);

    return 0;
}
