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

#include <igraph.h>
#include "test_utilities.h"

static igraph_int_t warning_count = 0;

static void counting_warning_handler(const char *reason, const char *file, int line) {
    (void) reason;
    (void) file;
    (void) line;
    warning_count++;
}

/* Compares the sparse result against the dense version, including the row and column ids. */
static void check_against_dense(
    const igraph_t *graph, const igraph_vector_bool_t *types,
    const igraph_vector_t *weights
) {
    igraph_matrix_t dense, from_sparse;
    igraph_sparsemat_t sparse;
    igraph_vector_int_t row_ids_dense, col_ids_dense, row_ids_sparse, col_ids_sparse;

    igraph_matrix_init(&dense, 0, 0);
    igraph_matrix_init(&from_sparse, 0, 0);
    igraph_sparsemat_init(&sparse, 0, 0, 0);
    igraph_vector_int_init(&row_ids_dense, 0);
    igraph_vector_int_init(&col_ids_dense, 0);
    igraph_vector_int_init(&row_ids_sparse, 0);
    igraph_vector_int_init(&col_ids_sparse, 0);

    IGRAPH_ASSERT(igraph_get_biadjacency(
        graph, types, weights, &dense, &row_ids_dense, &col_ids_dense) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_get_biadjacency_sparse(
        graph, types, weights, &sparse, &row_ids_sparse, &col_ids_sparse) == IGRAPH_SUCCESS);

    IGRAPH_ASSERT(igraph_sparsemat_as_matrix(&from_sparse, &sparse) == IGRAPH_SUCCESS);

    IGRAPH_ASSERT(igraph_matrix_all_e(&dense, &from_sparse));
    IGRAPH_ASSERT(igraph_vector_int_all_e(&row_ids_dense, &row_ids_sparse));
    IGRAPH_ASSERT(igraph_vector_int_all_e(&col_ids_dense, &col_ids_sparse));

    printf("Sparse biadjacency matrix:\n");
    print_matrix(&from_sparse);
    printf("Row ids:\n");
    print_vector_int(&row_ids_sparse);
    printf("Column ids:\n");
    print_vector_int(&col_ids_sparse);
    printf("\n");

    igraph_vector_int_destroy(&col_ids_sparse);
    igraph_vector_int_destroy(&row_ids_sparse);
    igraph_vector_int_destroy(&col_ids_dense);
    igraph_vector_int_destroy(&row_ids_dense);
    igraph_sparsemat_destroy(&sparse);
    igraph_matrix_destroy(&from_sparse);
    igraph_matrix_destroy(&dense);
}

int main(void) {
    igraph_t graph;
    igraph_vector_bool_t types;
    igraph_vector_t weights;

    igraph_set_warning_handler(counting_warning_handler);

    /* Bipartite graph on 6 vertices. Vertices 0, 2, 4 form the first partition, 1, 3, 5 the second.
     * The edges 0-1 appear twice, to check that edge multiplicities are added. */
    igraph_small(&graph, 6, IGRAPH_UNDIRECTED,
                 0, 1, 0, 1, 0, 3, 2, 3, 2, 5, 4, 1, 4, 5,
                 -1);

    igraph_vector_bool_init_int(&types, 6, 0, 1, 0, 1, 0, 1);
    igraph_vector_init(&weights, 7);
    for (igraph_int_t i = 0; i < 7; i++) {
        VECTOR(weights)[i] = 0.5 + (igraph_real_t) i;
    }

    printf("Unweighted:\n");
    check_against_dense(&graph, &types, NULL);
    IGRAPH_ASSERT(warning_count == 0);

    printf("Weighted:\n");
    check_against_dense(&graph, &types, &weights);
    IGRAPH_ASSERT(warning_count == 0);

    igraph_vector_destroy(&weights);
    igraph_vector_bool_destroy(&types);
    igraph_destroy(&graph);

    /* An edge within the first partition must be ignored and must trigger a warning from both the dense and the sparse function. */
    igraph_small(&graph, 4, IGRAPH_UNDIRECTED, 0, 2, 1, 3, 0, 1, -1);
    igraph_vector_bool_init_int(&types, 4, 0, 0, 1, 1);

    printf("Edge within a partition:\n");
    warning_count = 0;
    check_against_dense(&graph, &types, NULL);
    IGRAPH_ASSERT(warning_count == 2);

    /* Wrong type vector size must be rejected. */
    igraph_vector_bool_resize(&types, 3);
    igraph_sparsemat_t sparse;
    igraph_sparsemat_init(&sparse, 0, 0, 0);
    CHECK_ERROR(igraph_get_biadjacency_sparse(&graph, &types, NULL, &sparse, NULL, NULL), IGRAPH_EINVAL);
    igraph_sparsemat_destroy(&sparse);

    /* Wrong weight vector size must be rejected. */
    igraph_vector_bool_resize(&types, 4);
    igraph_vector_init_int(&weights, 2, 1, 1);
    igraph_sparsemat_init(&sparse, 0, 0, 0);
    CHECK_ERROR(igraph_get_biadjacency_sparse(&graph, &types, &weights, &sparse, NULL, NULL), IGRAPH_EINVAL);
    igraph_sparsemat_destroy(&sparse);
    igraph_vector_destroy(&weights);

    igraph_vector_bool_destroy(&types);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
