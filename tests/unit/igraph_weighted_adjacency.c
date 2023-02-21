/* -*- mode: C -*-  */
/*
    IGraph library.
    Copyright (C) 2006-2022  The igraph development team

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
#include <stdarg.h>

#include "test_utilities.h"

//TODO make public?
igraph_error_t igraph_is_same_graph_weighted(const igraph_t *graph1, const igraph_t *graph2, igraph_bool_t *res, igraph_vector_t *weights1, igraph_vector_t *weights2) {
    igraph_integer_t nv1 = igraph_vcount(graph1);
    igraph_integer_t nv2 = igraph_vcount(graph2);
    igraph_integer_t ne1 = igraph_ecount(graph1);
    igraph_integer_t ne2 = igraph_ecount(graph2);
    igraph_integer_t i, eid1, eid2;

    *res = 0; /* Assume that the graphs differ */

    /* Check for same number of vertices/edges */
    if ((nv1 != nv2) || (ne1 != ne2)) {
        return IGRAPH_SUCCESS;
    }

    /* Check for same directedness */
    if (igraph_is_directed(graph1) != igraph_is_directed(graph2)) {
        return IGRAPH_SUCCESS;
    }

    for (i = 0; i < ne1; i++) {
        eid1 = VECTOR(graph1->ii)[i];
        eid2 = VECTOR(graph2->ii)[i];

        /* Check they have the same source */
        if (IGRAPH_FROM(graph1, eid1) != IGRAPH_FROM(graph2, eid2)) {
            return IGRAPH_SUCCESS;
        }

        /* Check they have the same target */
        if (IGRAPH_TO(graph1, eid1) != IGRAPH_TO(graph2, eid2)) {
            return IGRAPH_SUCCESS;
        }

        /* Check they have the same weight */
        if (VECTOR(*weights1)[eid1] != VECTOR(*weights2)[eid2]) {
            return IGRAPH_SUCCESS;
        }
    }

    *res = 1; /* No difference was found, graphs are the same */
    return IGRAPH_SUCCESS;
}

void check_sparsemat(igraph_t *g, igraph_adjacency_t mode, igraph_loops_t loops, igraph_matrix_t *adjmatrix, igraph_vector_t *other_weights) {
    igraph_t g_sparse;
    igraph_sparsemat_t sparse_adjmatrix, sparse_adjmatrix_comp;
    igraph_bool_t same;
    igraph_vector_t weights;

    igraph_vector_init(&weights, 0);

    igraph_matrix_as_sparsemat(&sparse_adjmatrix, adjmatrix, 0.0001);
    igraph_sparsemat_compress(&sparse_adjmatrix, &sparse_adjmatrix_comp);
    igraph_sparse_weighted_adjacency(&g_sparse, &sparse_adjmatrix_comp, mode, &weights, loops);
    igraph_is_same_graph_weighted(g, &g_sparse, &same, other_weights, &weights);
    if (!same) {
        printf("Sparse graph differs from non-sparse:\n");
        print_weighted_graph(&g_sparse, &weights);
        exit(1);
    }
    igraph_sparsemat_destroy(&sparse_adjmatrix);
    igraph_sparsemat_destroy(&sparse_adjmatrix_comp);
    igraph_destroy(&g_sparse);
    igraph_vector_destroy(&weights);
}

void test_single_matrix(igraph_matrix_t* mat, igraph_adjacency_t mode) {
    igraph_t g;
    igraph_vector_t weights;

    igraph_vector_init(&weights, 0);

    printf("No loops\n--------\n\n");
    igraph_weighted_adjacency(&g, mat, mode, &weights, IGRAPH_NO_LOOPS);

    print_weighted_graph(&g, &weights);
    check_sparsemat(&g, mode, IGRAPH_NO_LOOPS, mat, &weights);

    igraph_destroy(&g);
    printf("\n");

    printf("Loops once\n----------\n\n");
    igraph_weighted_adjacency(&g, mat, mode, &weights, IGRAPH_LOOPS_ONCE);
    print_weighted_graph(&g, &weights);
    check_sparsemat(&g, mode, IGRAPH_LOOPS_ONCE, mat, &weights);
    igraph_destroy(&g);
    printf("\n");

    printf("Loops twice\n-----------\n\n");
    igraph_weighted_adjacency(&g, mat, mode, &weights, IGRAPH_LOOPS_TWICE);
    print_weighted_graph(&g, &weights);
    check_sparsemat(&g, mode, IGRAPH_LOOPS_TWICE, mat, &weights);
    igraph_destroy(&g);
    printf("\n");

    igraph_vector_destroy(&weights);

    VERIFY_FINALLY_STACK();
}

void check_error(igraph_matrix_t *adjmatrix, igraph_adjacency_t mode, igraph_loops_t loops, igraph_error_t error) {
    igraph_t g;
    igraph_sparsemat_t sparse_adjmatrix, sparse_adjmatrix_comp;
    igraph_vector_t weights;
    igraph_vector_init(&weights, 0);

    igraph_matrix_as_sparsemat(&sparse_adjmatrix, adjmatrix, 0.0001);
    igraph_sparsemat_compress(&sparse_adjmatrix, &sparse_adjmatrix_comp);

    CHECK_ERROR(igraph_weighted_adjacency(&g, adjmatrix, mode, &weights, loops), error);
    CHECK_ERROR(igraph_sparse_adjacency(&g, &sparse_adjmatrix_comp, mode, loops), error);

    igraph_sparsemat_destroy(&sparse_adjmatrix);
    igraph_sparsemat_destroy(&sparse_adjmatrix_comp);
    igraph_vector_destroy(&weights);
}

int main(void) {
    igraph_matrix_t mat;
    igraph_matrix_t mat_sym;
    int m[4][4] = { { 0, 1, 2, 0 }, { 2, 0, 0, 1 }, { 0, 0, 4, 0 }, { 0, 1, 0, 0 } };
    int m_sym[4][4] = { { 0, 2, 2, 0 }, { 2, 0, 0, 1 }, { 2, 0, 4, 0 }, { 0, 1, 0, 0 } };
    igraph_integer_t i, j;

    printf("0x0 matrix\n");
    printf("==========\n\n");
    igraph_matrix_init(&mat, 0, 0);
    test_single_matrix(&mat, IGRAPH_ADJ_DIRECTED);
    igraph_matrix_destroy(&mat);

    igraph_matrix_init(&mat, 4, 4);
    for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) {
        MATRIX(mat, i, j) = m[i][j];
    }
    igraph_matrix_init(&mat_sym, 4, 4);
    for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) {
        MATRIX(mat_sym, i, j) = m_sym[i][j];
    }

    /* [ 0 1 2 0 ]
       [ 2 0 0 1 ]
       [ 0 0 4 0 ]
       [ 0 1 0 0 ] */
    printf("IGRAPH_ADJ_DIRECTED\n");
    printf("===================\n\n");
    test_single_matrix(&mat, IGRAPH_ADJ_DIRECTED);

    /* [ 0 1 2 0 ]
       [ - 0 0 1 ]
       [ - - 4 0 ]
       [ - - - 0 ] */
    printf("IGRAPH_ADJ_UPPER\n");
    printf("================\n\n");
    test_single_matrix(&mat, IGRAPH_ADJ_UPPER);

    /* [ 0 - - - ]
       [ 2 0 - - ]
       [ 0 0 4 - ]
       [ 0 1 0 0 ] */
    printf("IGRAPH_ADJ_LOWER\n");
    printf("================\n\n");
    test_single_matrix(&mat, IGRAPH_ADJ_LOWER);

    /* [ 0 1 0 0 ]
       [ 1 0 0 1 ]
       [ 0 0 4 0 ]
       [ 0 1 0 0 ] */
    printf("IGRAPH_ADJ_MIN\n");
    printf("==============\n\n");
    test_single_matrix(&mat, IGRAPH_ADJ_MIN);

    /* [ 0 2 2 0 ]
       [ 2 0 0 1 ]
       [ 2 0 4 0 ]
       [ 0 1 0 0 ] */
    printf("IGRAPH_ADJ_MAX\n");
    printf("==============\n\n");
    test_single_matrix(&mat, IGRAPH_ADJ_MAX);

    /* [ 0 2 2 0 ]
       [ 2 0 0 1 ]
       [ 2 0 4 0 ]
       [ 0 1 0 0 ] */
    printf("IGRAPH_ADJ_UNDIRECTED\n");
    printf("==============\n\n");
    test_single_matrix(&mat_sym, IGRAPH_ADJ_UNDIRECTED);

    /* [ 0 3 2 0 ]
       [ 3 0 0 2 ]
       [ 2 0 4 0 ]
       [ 0 2 0 0 ] */
    printf("IGRAPH_ADJ_PLUS\n");
    printf("===============\n\n");
    test_single_matrix(&mat, IGRAPH_ADJ_PLUS);

    igraph_matrix_destroy(&mat);
    igraph_matrix_destroy(&mat_sym);

    VERIFY_FINALLY_STACK();

    {
        printf("Check handling of non-square matrix error.\n");
        igraph_real_t e[] = {1, 2, 0};
        igraph_matrix_view(&mat, e, 1, 3);
        check_error(&mat, IGRAPH_ADJ_DIRECTED, IGRAPH_NO_LOOPS, IGRAPH_NONSQUARE);
    }

    {
        printf("Check handling of invalid adjacency mode.\n");
        igraph_real_t e[] = {0, 2, 0, 3, 0, 4, 0, 5, 6};
        igraph_matrix_view(&mat, e, 3, 3);
        check_error(&mat, (igraph_adjacency_t) 42, IGRAPH_LOOPS_TWICE, IGRAPH_EINVAL);
    }

    {
        printf("Check error for 0x1 matrix.\n");
        igraph_matrix_init(&mat, 0, 1);
        check_error(&mat, IGRAPH_ADJ_DIRECTED, IGRAPH_LOOPS_TWICE, IGRAPH_NONSQUARE);
        igraph_matrix_destroy(&mat);
    }

    {
        printf("Check error for non-symmetric matrix and IGRAPH_ADJ_UNDIRECTED.\n");
        igraph_real_t e[] = {0, 2, 0, 3, 0, 4, 0, 5, 6};
        igraph_matrix_view(&mat, e, 3, 3);
        check_error(&mat, IGRAPH_ADJ_UNDIRECTED, IGRAPH_LOOPS_TWICE, IGRAPH_EINVAL);
    }

    VERIFY_FINALLY_STACK();

    return 0;
}
