/*
   igraph library.
   Copyright (C) 2005-2021 The igraph development team

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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "igraph_constructors.h"

#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_sparsemat.h"

static igraph_error_t igraph_i_adjacency_directed_or_plus(
    const igraph_matrix_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_adjacency_t mode, igraph_loops_t loops
);
static igraph_error_t igraph_i_adjacency_max(
    const igraph_matrix_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_loops_t loops
);
static igraph_error_t igraph_i_adjacency_undirected(
    const igraph_matrix_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_loops_t loops
);
static igraph_error_t igraph_i_adjacency_upper(
    const igraph_matrix_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_loops_t loops
);
static igraph_error_t igraph_i_adjacency_lower(
    const igraph_matrix_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_loops_t loops
);
static igraph_error_t igraph_i_adjacency_min(
    const igraph_matrix_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_loops_t loops
);

static igraph_error_t igraph_i_adjust_loop_edge_count(
    igraph_int_t* count, igraph_loops_t loops
) {
    /* The compiler should be smart enough to figure out that this can be
     * inlined */
    switch (loops) {
        case IGRAPH_NO_LOOPS:
            *count = 0;
            break;
        case IGRAPH_LOOPS_TWICE:
            if (*count & 1) {
                IGRAPH_ERROR("Odd number found in the diagonal of the adjacency matrix.", IGRAPH_EINVAL);
            }
            *count >>= 1;
            break;
        case IGRAPH_LOOPS_ONCE:
        default:
            break;
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_adjacency_directed_or_plus(
    const igraph_matrix_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_adjacency_t mode, igraph_loops_t loops
) {
    const igraph_int_t no_of_nodes = igraph_matrix_nrow(adjmatrix);

    /* For sake of consistency with the rest of the library, IGRAPH_LOOPS_TWICE
     * is treated as IGRAPH_LOOPS_ONCE for directed graphs */
    if (mode == IGRAPH_ADJ_DIRECTED && loops == IGRAPH_LOOPS_TWICE) {
        loops = IGRAPH_LOOPS_ONCE;
    }

    for (igraph_int_t j = 0; j < no_of_nodes; j++) {
        for (igraph_int_t i = 0; i < no_of_nodes; i++) {
            igraph_int_t M = MATRIX(*adjmatrix, i, j);
            if (i == j) {
                IGRAPH_CHECK(igraph_i_adjust_loop_edge_count(&M, loops));
            }
            for (igraph_int_t k = 0; k < M; k++) {
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, j));
            }
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_adjacency_max(
    const igraph_matrix_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_loops_t loops
) {

    igraph_int_t no_of_nodes = igraph_matrix_nrow(adjmatrix);
    igraph_int_t i, j, k, M1, M2;

    for (i = 0; i < no_of_nodes; i++) {
        /* do the loops first */
        M1 = MATRIX(*adjmatrix, i, i);
        if (M1) {
            IGRAPH_CHECK(igraph_i_adjust_loop_edge_count(&M1, loops));
            for (k = 0; k < M1; k++) {
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, i));
            }
        }

        /* then the rest */
        for (j = i + 1; j < no_of_nodes; j++) {
            M1 = MATRIX(*adjmatrix, i, j);
            M2 = MATRIX(*adjmatrix, j, i);
            if (M1 < M2) {
                M1 = M2;
            }
            for (k = 0; k < M1; k++) {
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, j));
            }
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_adjacency_undirected(
    const igraph_matrix_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_loops_t loops
) {
    if (!igraph_matrix_is_symmetric(adjmatrix)) {
        IGRAPH_ERROR(
            "Adjacency matrix should be symmetric to produce an undirected graph.",
            IGRAPH_EINVAL
        );
    }
    return igraph_i_adjacency_max(adjmatrix, edges, loops);
}

static igraph_error_t igraph_i_adjacency_upper(
    const igraph_matrix_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_loops_t loops
) {

    const igraph_int_t no_of_nodes = igraph_matrix_nrow(adjmatrix);
    igraph_int_t M;

    /* IGRAPH_LOOPS_TWICE is treated as IGRAPH_LOOPS_ONCE -- it makes no sense
     * for loops to appear twice in the adjacency matrix when the lower triangle
     * is empty; double-counting of loops in undirected graphs happens because
     * the upper and the lower triangle are added on top of each other on the
     * diagonal. See discussion in #2501:
     *
     * https://github.com/igraph/igraph/issues/2501#issuecomment-1949345675 */
    if (loops == IGRAPH_LOOPS_TWICE) {
        loops = IGRAPH_LOOPS_ONCE;
    }

    for (igraph_int_t j = 0; j < no_of_nodes; j++) {
        for (igraph_int_t i = 0; i < j; i++) {
            M = MATRIX(*adjmatrix, i, j);
            for (igraph_int_t k = 0; k < M; k++) {
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, j));
            }
        }

        /* do the loops as well */
        M = MATRIX(*adjmatrix, j, j);
        if (M) {
            IGRAPH_CHECK(igraph_i_adjust_loop_edge_count(&M, loops));
            for (igraph_int_t k = 0; k < M; k++) {
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, j));
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, j));
            }
        }
    }
    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_adjacency_lower(
    const igraph_matrix_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_loops_t loops
) {

    const igraph_int_t no_of_nodes = igraph_matrix_nrow(adjmatrix);
    igraph_int_t M;

    /* IGRAPH_LOOPS_TWICE is treated as IGRAPH_LOOPS_ONCE -- it makes no sense
     * for loops to appear twice in the adjacency matrix when the lower triangle
     * is empty; double-counting of loops in undirected graphs happens because
     * the upper and the lower triangle are added on top of each other on the
     * diagonal. See discussion in #2501:
     *
     * https://github.com/igraph/igraph/issues/2501#issuecomment-1949345675 */
    if (loops == IGRAPH_LOOPS_TWICE) {
        loops = IGRAPH_LOOPS_ONCE;
    }

    for (igraph_int_t j = 0; j < no_of_nodes; j++) {
        /* do the loops first */
        M = MATRIX(*adjmatrix, j, j);
        if (M) {
            IGRAPH_CHECK(igraph_i_adjust_loop_edge_count(&M, loops));
            for (igraph_int_t k = 0; k < M; k++) {
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, j));
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, j));
            }
        }

        for (igraph_int_t i = j+1; i < no_of_nodes; i++) {
            M = MATRIX(*adjmatrix, i, j);
            for (igraph_int_t k = 0; k < M; k++) {
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, j));
            }
        }
    }
    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_adjacency_min(
    const igraph_matrix_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_loops_t loops
) {

    igraph_int_t no_of_nodes = igraph_matrix_nrow(adjmatrix);
    igraph_int_t i, j, k, M1, M2;

    for (i = 0; i < no_of_nodes; i++) {
        /* do the loops first */
        M1 = MATRIX(*adjmatrix, i, i);
        if (M1) {
            IGRAPH_CHECK(igraph_i_adjust_loop_edge_count(&M1, loops));
            for (k = 0; k < M1; k++) {
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, i));
            }
        }

        /* then the rest */
        for (j = i + 1; j < no_of_nodes; j++) {
            M1 = MATRIX(*adjmatrix, i, j);
            M2 = MATRIX(*adjmatrix, j, i);
            if (M1 > M2) {
                M1 = M2;
            }
            for (k = 0; k < M1; k++) {
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, j));
            }
        }
    }

    return IGRAPH_SUCCESS;
}



/**
 * \ingroup generators
 * \function igraph_adjacency
 * \brief Creates a graph from an adjacency matrix.
 *
 * The order of the vertices in the matrix is preserved, i.e. the vertex
 * corresponding to the first row/column will be vertex with id 0, the
 * next row is for vertex 1, etc. No guarantees are given about the ordering
 * of edges.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param adjmatrix The adjacency matrix. How it is interpreted
 *        depends on the \p mode argument.
 * \param mode Constant to specify how the given matrix is interpreted
 *        as an adjacency matrix. Possible values (A(i,j) is the element in
 *        row i and column j in the adjacency matrix \p adjmatrix):
 *        \clist
 *        \cli IGRAPH_ADJ_DIRECTED
 *          The graph will be directed and an element gives the number of edges
 *           between two vertices.
 *        \cli IGRAPH_ADJ_UNDIRECTED
 *          The graph will be undirected and an element gives the number of
 *          edges between two vertices. If the input matrix is not symmetric,
 *          an error is thrown.
 *        \cli IGRAPH_ADJ_MAX
 *          An undirected graph will be created and the number of edges between
 *          vertices i and j is max(A(i,j), A(j,i)).
 *        \cli IGRAPH_ADJ_MIN
 *          An undirected graph will be created with min(A(i,j), A(j,i)) edges
 *          between vertices i and j.
 *        \cli IGRAPH_ADJ_PLUS
 *          An undirected graph will be created with A(i,j)+A(j,i) edges
 *          between vertices i and j.
 *        \cli IGRAPH_ADJ_UPPER
 *          An undirected graph will be created. Only the upper right triangle
 *          (including the diagonal) is used for the number of edges.
 *        \cli IGRAPH_ADJ_LOWER
 *          An undirected graph will be created. Only the lower left triangle
 *          (including the diagonal) is used for the number of edges.
 *       \endclist
 * \param loops Constant of type \ref igraph_loops_t to specify how the diagonal
 *        of the matrix should be treated when creating loop edges. Ignored for
 *        modes \c IGRAPH_ADJ_DIRECTED, \c IGRAPH_ADJ_UPPER and \c IGRAPH_ADJ_LOWER.
 *        \clist
 *        \cli IGRAPH_NO_LOOPS
 *          Ignore the diagonal of the input matrix and do not create loops.
 *        \cli IGRAPH_LOOPS_ONCE
 *          Treat the diagonal entries as the number of loop edges incident on
 *          the corresponding vertex.
 *        \cli IGRAPH_LOOPS_TWICE
 *          Treat the diagonal entries as \em twice the number of loop edges
 *          incident on the corresponding vertex. Odd numbers in the diagonal
 *          will return an error code.
 *        \endclist
 * \return Error code,
 *         \c IGRAPH_EINVAL: Non-square adjacency matrix, negative entry in
 *         adjacency matrix, or an odd number was found in the diagonal with
 *         \c IGRAPH_LOOPS_TWICE
 *
 * Time complexity: O(|V||V|),
 * |V| is the number of vertices in the graph.
 */
igraph_error_t igraph_adjacency(
    igraph_t *graph, const igraph_matrix_t *adjmatrix, igraph_adjacency_t mode,
    igraph_loops_t loops
) {

    igraph_vector_int_t edges = IGRAPH_VECTOR_NULL;
    igraph_int_t no_of_nodes = igraph_matrix_nrow(adjmatrix);

    /* Some checks */
    if (igraph_matrix_nrow(adjmatrix) != igraph_matrix_ncol(adjmatrix)) {
        IGRAPH_ERROR("Adjacency matrices must be square.", IGRAPH_EINVAL);
    }

    if (no_of_nodes != 0 && igraph_matrix_min(adjmatrix) < 0) {
        IGRAPH_ERRORF("Edge counts should be non-negative, found %g.", IGRAPH_EINVAL,
                igraph_matrix_min(adjmatrix));
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

    /* Collect the edges */
    no_of_nodes = igraph_matrix_nrow(adjmatrix);
    switch (mode) {
    case IGRAPH_ADJ_DIRECTED:
    case IGRAPH_ADJ_PLUS:
        IGRAPH_CHECK(igraph_i_adjacency_directed_or_plus(adjmatrix, &edges, mode, loops));
        break;
    case IGRAPH_ADJ_MAX:
        IGRAPH_CHECK(igraph_i_adjacency_max(adjmatrix, &edges, loops));
        break;
    case IGRAPH_ADJ_UNDIRECTED:
        IGRAPH_CHECK(igraph_i_adjacency_undirected(adjmatrix, &edges, loops));
        break;
    case IGRAPH_ADJ_UPPER:
        IGRAPH_CHECK(igraph_i_adjacency_upper(adjmatrix, &edges, loops));
        break;
    case IGRAPH_ADJ_LOWER:
        IGRAPH_CHECK(igraph_i_adjacency_lower(adjmatrix, &edges, loops));
        break;
    case IGRAPH_ADJ_MIN:
        IGRAPH_CHECK(igraph_i_adjacency_min(adjmatrix, &edges, loops));
        break;
    default:
        IGRAPH_ERROR("Invalid adjacency mode.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, (mode == IGRAPH_ADJ_DIRECTED)));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_weighted_adjacency_directed(
        const igraph_matrix_t *adjmatrix,
        igraph_vector_int_t *edges,
        igraph_vector_t *weights,
        igraph_loops_t loops
);
static igraph_error_t igraph_i_weighted_adjacency_plus(
        const igraph_matrix_t *adjmatrix,
        igraph_vector_int_t *edges,
        igraph_vector_t *weights,
        igraph_loops_t loops
);
static igraph_error_t igraph_i_weighted_adjacency_max(
        const igraph_matrix_t *adjmatrix,
        igraph_vector_int_t *edges,
        igraph_vector_t *weights,
        igraph_loops_t loops
);
static igraph_error_t igraph_i_weighted_adjacency_upper(
        const igraph_matrix_t *adjmatrix,
        igraph_vector_int_t *edges,
        igraph_vector_t *weights,
        igraph_loops_t loops
);
static igraph_error_t igraph_i_weighted_adjacency_lower(
        const igraph_matrix_t *adjmatrix,
        igraph_vector_int_t *edges,
        igraph_vector_t *weights,
        igraph_loops_t loops
);
static igraph_error_t igraph_i_weighted_adjacency_min(
        const igraph_matrix_t *adjmatrix,
        igraph_vector_int_t *edges,
        igraph_vector_t *weights,
        igraph_loops_t loops
);

static void igraph_i_adjust_loop_edge_weight(igraph_real_t* weight, igraph_loops_t loops) {
    /* The compiler should be smart enough to figure out that this can be
     * inlined */
    switch (loops) {
        case IGRAPH_NO_LOOPS:
            *weight = 0.0;
            break;
        case IGRAPH_LOOPS_TWICE:
            *weight /= 2;
            break;
        case IGRAPH_LOOPS_ONCE:
        default:
            break;
    }
}

static igraph_error_t igraph_i_weighted_adjacency_directed(
    const igraph_matrix_t *adjmatrix,
    igraph_vector_int_t *edges,
    igraph_vector_t *weights,
    igraph_loops_t loops
) {

    const igraph_int_t no_of_nodes = igraph_matrix_nrow(adjmatrix);

    /* For sake of consistency with the rest of the library, IGRAPH_LOOPS_TWICE
     * is treated as IGRAPH_LOOPS_ONCE for directed graphs */
    if (loops == IGRAPH_LOOPS_TWICE) {
        loops = IGRAPH_LOOPS_ONCE;
    }

    for (igraph_int_t j = 0; j < no_of_nodes; j++) {
        for (igraph_int_t i = 0; i < no_of_nodes; i++) {
            igraph_real_t M = MATRIX(*adjmatrix, i, j);
            if (M != 0.0) {
                if (i == j) {
                    igraph_i_adjust_loop_edge_weight(&M, loops);
                    if (M == 0.0) {
                        continue;
                    }
                }
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, j));
                IGRAPH_CHECK(igraph_vector_push_back(weights, M));
            }
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_weighted_adjacency_plus(
    const igraph_matrix_t *adjmatrix,
    igraph_vector_int_t *edges,
    igraph_vector_t *weights,
    igraph_loops_t loops
) {

    igraph_int_t no_of_nodes = igraph_matrix_nrow(adjmatrix);
    igraph_int_t i, j;
    igraph_real_t M;

    for (i = 0; i < no_of_nodes; i++) {
        if (loops != IGRAPH_NO_LOOPS) {
            M = MATRIX(*adjmatrix, i, i);
            if (M != 0.0) {
                igraph_i_adjust_loop_edge_weight(&M, loops);
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(weights, M));
            }
        }

        for (j = i + 1; j < no_of_nodes; j++) {
            M = MATRIX(*adjmatrix, i, j) + MATRIX(*adjmatrix, j, i);
            if (M != 0.0) {
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, j));
                IGRAPH_CHECK(igraph_vector_push_back(weights, M));
            }
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_weighted_adjacency_max(
    const igraph_matrix_t *adjmatrix,
    igraph_vector_int_t *edges,
    igraph_vector_t *weights,
    igraph_loops_t loops
) {

    igraph_int_t no_of_nodes = igraph_matrix_nrow(adjmatrix);
    igraph_int_t i, j;
    igraph_real_t M1, M2;

    for (i = 0; i < no_of_nodes; i++) {
        /* do the loops first */
        if (loops) {
            M1 = MATRIX(*adjmatrix, i, i);
            if (M1 != 0.0) {
                igraph_i_adjust_loop_edge_weight(&M1, loops);
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(weights, M1));
            }
        }

        /* then the rest */
        for (j = i + 1; j < no_of_nodes; j++) {
            M1 = MATRIX(*adjmatrix, i, j);
            M2 = MATRIX(*adjmatrix, j, i);
            if (M1 < M2 || isnan(M2)) {
                M1 = M2;
            }
            if (M1 != 0.0) {
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, j));
                IGRAPH_CHECK(igraph_vector_push_back(weights, M1));
            }
        }
    }
    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_weighted_adjacency_undirected(
    const igraph_matrix_t *adjmatrix,
    igraph_vector_int_t *edges,
    igraph_vector_t *weights,
    igraph_loops_t loops
) {
    /* We do not use igraph_matrix_is_symmetric() for this check, as we need to
     * allow symmetric matrices with NaN values. igraph_matrix_is_symmetric()
     * returns false for these as NaN != NaN. */
    const igraph_int_t n = igraph_matrix_nrow(adjmatrix);
    for (igraph_int_t i=0; i < n; i++) {
        /* do the loops first */
        if (loops) {
            igraph_real_t M = MATRIX(*adjmatrix, i, i);
            if (M != 0.0) {
                igraph_i_adjust_loop_edge_weight(&M, loops);
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(weights, M));
            }
        }

        for (igraph_int_t j=0; j < i; j++) {
            igraph_real_t M1 = MATRIX(*adjmatrix, i, j);
            igraph_real_t M2 = MATRIX(*adjmatrix, j, i);
            if (IGRAPH_UNLIKELY(M1 != M2 && ! (isnan(M1) && isnan(M2)))) {
                IGRAPH_ERROR(
                    "Adjacency matrix should be symmetric to produce an undirected graph.",
                    IGRAPH_EINVAL
                );
            } else if (M1 != 0.0) {
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, j));
                IGRAPH_CHECK(igraph_vector_push_back(weights, M1));
            }
        }
    }
    return IGRAPH_SUCCESS;
}


static igraph_error_t igraph_i_weighted_adjacency_upper(
    const igraph_matrix_t *adjmatrix,
    igraph_vector_int_t *edges,
    igraph_vector_t *weights,
    igraph_loops_t loops
) {

    const igraph_int_t no_of_nodes = igraph_matrix_nrow(adjmatrix);
    igraph_real_t M;

    /* IGRAPH_LOOPS_TWICE is treated as IGRAPH_LOOPS_ONCE -- it makes no sense
     * for loops to appear twice in the adjacency matrix when the lower triangle
     * is empty; double-counting of loops in undirected graphs happens because
     * the upper and the lower triangle are added on top of each other on the
     * diagonal. See discussion in #2501:
     *
     * https://github.com/igraph/igraph/issues/2501#issuecomment-1949345675 */
    if (loops == IGRAPH_LOOPS_TWICE) {
        loops = IGRAPH_LOOPS_ONCE;
    }

    for (igraph_int_t j = 0; j < no_of_nodes; j++) {
        for (igraph_int_t i = 0; i < j; i++) {
            igraph_real_t M = MATRIX(*adjmatrix, i, j);
            if (M != 0.0) {
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, j));
                IGRAPH_CHECK(igraph_vector_push_back(weights, M));
            }
        }

        /* do the loops as well */
        if (loops) {
            M = MATRIX(*adjmatrix, j, j);
            if (M != 0.0) {
                igraph_i_adjust_loop_edge_weight(&M, loops);
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, j));
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, j));
                IGRAPH_CHECK(igraph_vector_push_back(weights, M));
            }
        }
    }
    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_weighted_adjacency_lower(
    const igraph_matrix_t *adjmatrix,
    igraph_vector_int_t *edges,
    igraph_vector_t *weights,
    igraph_loops_t loops
) {

    const igraph_int_t no_of_nodes = igraph_matrix_nrow(adjmatrix);
    igraph_real_t M;

    /* IGRAPH_LOOPS_TWICE is treated as IGRAPH_LOOPS_ONCE -- it makes no sense
     * for loops to appear twice in the adjacency matrix when the lower triangle
     * is empty; double-counting of loops in undirected graphs happens because
     * the upper and the lower triangle are added on top of each other on the
     * diagonal. See discussion in #2501:
     *
     * https://github.com/igraph/igraph/issues/2501#issuecomment-1949345675 */
    if (loops == IGRAPH_LOOPS_TWICE) {
        loops = IGRAPH_LOOPS_ONCE;
    }

    for (igraph_int_t j = 0; j < no_of_nodes; j++) {
        /* do the loops first */
        if (loops) {
            M = MATRIX(*adjmatrix, j, j);
            if (M != 0.0) {
                igraph_i_adjust_loop_edge_weight(&M, loops);
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, j));
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, j));
                IGRAPH_CHECK(igraph_vector_push_back(weights, M));
            }
        }

        for (igraph_int_t i = j+1; i < no_of_nodes; i++) {
            M = MATRIX(*adjmatrix, i, j);
            if (M != 0.0) {
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, j));
                IGRAPH_CHECK(igraph_vector_push_back(weights, M));
            }
        }
    }
    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_weighted_adjacency_min(
    const igraph_matrix_t *adjmatrix,
    igraph_vector_int_t *edges,
    igraph_vector_t *weights,
    igraph_loops_t loops
) {

    igraph_int_t no_of_nodes = igraph_matrix_nrow(adjmatrix);
    igraph_int_t i, j;
    igraph_real_t M1, M2;

    for (i = 0; i < no_of_nodes; i++) {
        /* do the loops first */
        if (loops) {
            M1 = MATRIX(*adjmatrix, i, i);
            if (M1 != 0.0) {
                igraph_i_adjust_loop_edge_weight(&M1, loops);
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(weights, M1));
            }
        }

        /* then the rest */
        for (j = i + 1; j < no_of_nodes; j++) {
            M1 = MATRIX(*adjmatrix, i, j);
            M2 = MATRIX(*adjmatrix, j, i);
            if (M1 > M2 || isnan(M2)) {
                M1 = M2;
            }
            if (M1 != 0.0) {
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, j));
                IGRAPH_CHECK(igraph_vector_push_back(weights, M1));
            }
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup generators
 * \function igraph_weighted_adjacency
 * \brief Creates a graph from a weighted adjacency matrix.
 *
 * The order of the vertices in the matrix is preserved, i.e. the vertex
 * corresponding to the first row/column will be vertex with id 0, the
 * next row is for vertex 1, etc. No guarantees are given for the ordering
 * of edges.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param adjmatrix The weighted adjacency matrix. How it is interpreted
 *        depends on the \p mode argument. The common feature is that
 *        edges with zero weights are considered nonexistent (however,
 *        negative weights are permitted).
 * \param mode Constant to specify how the given matrix is interpreted
 *        as an adjacency matrix. Possible values (A(i,j) is the element in row
 *        i and column j in the adjacency matrix \p adjmatrix):
 *        \clist
 *        \cli IGRAPH_ADJ_DIRECTED
 *          The graph will be directed and an element specifies the weight of the
 *           edge between two vertices.
 *        \cli IGRAPH_ADJ_UNDIRECTED
 *          This is the same as \c IGRAPH_ADJ_MAX, for convenience.
 *        \cli IGRAPH_ADJ_MAX
 *          An undirected graph will be created and the weight of the edge between
 *          vertices i and j is max(A(i,j), A(j,i)).
 *        \cli IGRAPH_ADJ_MIN
 *          An undirected graph will be created and the weight of the edge between
 *          vertices i and j is min(A(i,j), A(j,i)).
 *        \cli IGRAPH_ADJ_PLUS
 *          An undirected graph will be created and the weight of the edge between
 *          vertices i and j is A(i,j)+A(j,i).
 *        \cli IGRAPH_ADJ_UPPER
 *          An undirected graph will be created. Only the upper right triangle
 *          (including the diagonal) is used for the edge weights.
 *        \cli IGRAPH_ADJ_LOWER
 *          An undirected graph will be created. Only the lower left triangle
 *          (including the diagonal) is used for the edge weights.
 *       \endclist
 * \param weights Pointer to an initialized vector, the weights will be stored here.
 * \param loops Constant to specify how the diagonal of the matrix should be
 *        treated when creating loop edges. Ignored for modes
 *        \c IGRAPH_ADJ_DIRECTED, \c IGRAPH_ADJ_UPPER and \c IGRAPH_ADJ_LOWER.
 *        \clist
 *        \cli IGRAPH_NO_LOOPS
 *          Ignore the diagonal of the input matrix and do not create loops.
 *        \cli IGRAPH_LOOPS_ONCE
 *          Treat the diagonal entries as the weight of the loop edge incident
 *          on the corresponding vertex.
 *        \cli IGRAPH_LOOPS_TWICE
 *          Treat the diagonal entries as \em twice the weight of the loop edge
 *          incident on the corresponding vertex.
 *        \endclist
 * \return Error code,
 *         \c IGRAPH_EINVAL: non-square matrix.
 *
 * Time complexity: O(|V||V|),
 * |V| is the number of vertices in the graph.
 *
 * \example examples/simple/igraph_weighted_adjacency.c
 */
igraph_error_t igraph_weighted_adjacency(
    igraph_t *graph, const igraph_matrix_t *adjmatrix, igraph_adjacency_t mode,
    igraph_vector_t *weights, igraph_loops_t loops
) {

    igraph_vector_int_t edges = IGRAPH_VECTOR_NULL;
    igraph_int_t no_of_nodes;

    /* Some checks */
    if (igraph_matrix_nrow(adjmatrix) != igraph_matrix_ncol(adjmatrix)) {
        IGRAPH_ERROR("Adjacency matrices must be square.", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    igraph_vector_clear(weights);

    /* Collect the edges */
    no_of_nodes = igraph_matrix_nrow(adjmatrix);
    switch (mode) {
    case IGRAPH_ADJ_DIRECTED:
        IGRAPH_CHECK(igraph_i_weighted_adjacency_directed(adjmatrix, &edges,
                     weights, loops));
        break;
    case IGRAPH_ADJ_MAX:
        IGRAPH_CHECK(igraph_i_weighted_adjacency_max(adjmatrix, &edges,
                     weights, loops));
        break;
    case IGRAPH_ADJ_UNDIRECTED:
        IGRAPH_CHECK(igraph_i_weighted_adjacency_undirected(adjmatrix, &edges,
                     weights, loops));
        break;
    case IGRAPH_ADJ_UPPER:
        IGRAPH_CHECK(igraph_i_weighted_adjacency_upper(adjmatrix, &edges,
                     weights, loops));
        break;
    case IGRAPH_ADJ_LOWER:
        IGRAPH_CHECK(igraph_i_weighted_adjacency_lower(adjmatrix, &edges,
                     weights, loops));
        break;
    case IGRAPH_ADJ_MIN:
        IGRAPH_CHECK(igraph_i_weighted_adjacency_min(adjmatrix, &edges,
                     weights, loops));
        break;
    case IGRAPH_ADJ_PLUS:
        IGRAPH_CHECK(igraph_i_weighted_adjacency_plus(adjmatrix, &edges,
                     weights, loops));
        break;
    default:
        IGRAPH_ERROR("Invalid adjacency mode.", IGRAPH_EINVAL);
    }

    /* Create graph */
    IGRAPH_CHECK(igraph_empty(graph, no_of_nodes, (mode == IGRAPH_ADJ_DIRECTED)));
    IGRAPH_FINALLY(igraph_destroy, graph);
    if (igraph_vector_int_size(&edges) > 0) {
        IGRAPH_CHECK(igraph_add_edges(graph, &edges, NULL));
    }
    IGRAPH_FINALLY_CLEAN(1);

    /* Cleanup */
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_adjlist
 * \brief Creates a graph from an adjacency list.
 *
 * An adjacency list is a list of vectors, containing the neighbors
 * of all vertices. For operations that involve many changes to the
 * graph structure, it is recommended that you convert the graph into
 * an adjacency list via \ref igraph_adjlist_init(), perform the
 * modifications (these are cheap for an adjacency list) and then
 * recreate the igraph graph via this function.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param adjlist The adjacency list.
 * \param mode Whether or not to create a directed graph. \c IGRAPH_ALL
 *             means an undirected graph, \c IGRAPH_OUT means a
 *             directed graph from an out-adjacency list (i.e. each
 *             list contains the successors of the corresponding
 *             vertices), \c IGRAPH_IN means a directed graph from an
 *             in-adjacency list
 * \param duplicate Boolean constant. For undirected graphs this specifies
 *        whether each edge is included twice, in the vectors of
 *        both adjacent vertices. If this is \c false, then it is
 *        assumed that every edge is included only once. This argument
 *        is ignored for directed graphs.
 * \return Error code.
 *
 * \sa \ref igraph_adjlist_init() for the opposite operation.
 *
 * Time complexity: O(|V|+|E|).
 *
 */
igraph_error_t igraph_adjlist(igraph_t *graph, const igraph_adjlist_t *adjlist,
                   igraph_neimode_t mode, igraph_bool_t duplicate) {

    const igraph_int_t no_of_nodes = igraph_adjlist_size(adjlist);
    igraph_int_t no_of_edges = 0;

    igraph_vector_int_t edges;
    igraph_int_t edgeptr = 0;

    duplicate = duplicate && (mode == IGRAPH_ALL); /* only duplicate if undirected */

    for (igraph_int_t i = 0; i < no_of_nodes; i++) {
        no_of_edges += igraph_vector_int_size(igraph_adjlist_get(adjlist, i));
    }

    if (duplicate) {
        no_of_edges /= 2;
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 2 * no_of_edges);

    for (igraph_int_t i = 0; i < no_of_nodes; i++) {
        igraph_vector_int_t *neis = igraph_adjlist_get(adjlist, i);
        const igraph_int_t n = igraph_vector_int_size(neis);
        igraph_int_t loops = 0;

        for (igraph_int_t j = 0; j < n; j++) {
            igraph_int_t nei = VECTOR(*neis)[j];
            if (nei == i) {
                loops++;
            } else {
                if (! duplicate || nei > i) {
                    if (edgeptr + 2 > 2 * no_of_edges) {
                        IGRAPH_ERROR("Invalid adjacency list, most probably not correctly"
                                     " duplicated edges for an undirected graph.", IGRAPH_EINVAL);
                    }
                    if (mode == IGRAPH_IN) {
                        VECTOR(edges)[edgeptr++] = nei;
                        VECTOR(edges)[edgeptr++] = i;
                    } else {
                        VECTOR(edges)[edgeptr++] = i;
                        VECTOR(edges)[edgeptr++] = nei;
                    }
                }
            }
        }
        /* loops */
        if (duplicate) {
            loops = loops / 2;
        }
        if (edgeptr + 2 * loops > 2 * no_of_edges) {
            IGRAPH_ERROR("Invalid adjacency list, most probably not correctly"
                         " duplicated edges for an undirected graph.", IGRAPH_EINVAL);
        }
        for (igraph_int_t j = 0; j < loops; j++) {
            VECTOR(edges)[edgeptr++] = i;
            VECTOR(edges)[edgeptr++] = i;
        }
    }

    if (mode == IGRAPH_ALL) {
        IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, IGRAPH_UNDIRECTED));
    } else {
        IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, IGRAPH_DIRECTED));
    }

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_sparse_adjacency_directed_or_plus(
    igraph_sparsemat_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_adjacency_t mode, igraph_loops_t loops
) {
    igraph_sparsemat_iterator_t it;
    igraph_sparsemat_iterator_init(&it, adjmatrix);

    /* For sake of consistency with the rest of the library, IGRAPH_LOOPS_TWICE
     * is treated as IGRAPH_LOOPS_ONCE for directed graphs */
    if (mode == IGRAPH_ADJ_DIRECTED && loops == IGRAPH_LOOPS_TWICE) {
        loops = IGRAPH_LOOPS_ONCE;
    }

    for (; !igraph_sparsemat_iterator_end(&it); igraph_sparsemat_iterator_next(&it)) {
        igraph_int_t from = igraph_sparsemat_iterator_row(&it);
        igraph_int_t to = igraph_sparsemat_iterator_col(&it);
        igraph_int_t multi = igraph_sparsemat_iterator_get(&it);
        if (to == from) {
            IGRAPH_CHECK(igraph_i_adjust_loop_edge_count(&multi, loops));
        }
        for (igraph_int_t count = 0; count < multi; count++) {
            IGRAPH_CHECK(igraph_vector_int_push_back(edges, from));
            IGRAPH_CHECK(igraph_vector_int_push_back(edges, to));
        }
    }
    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_sparse_adjacency_max(
    igraph_sparsemat_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_loops_t loops
) {
    igraph_sparsemat_iterator_t it;
    igraph_real_t other;

    igraph_sparsemat_iterator_init(&it, adjmatrix);
    for (; !igraph_sparsemat_iterator_end(&it); igraph_sparsemat_iterator_next(&it)) {
        igraph_int_t from = igraph_sparsemat_iterator_row(&it);
        igraph_int_t to = igraph_sparsemat_iterator_col(&it);
        if (to < from) {
            continue;
        }
        igraph_int_t multi = igraph_sparsemat_iterator_get(&it);
        if (to == from) {
            IGRAPH_CHECK(igraph_i_adjust_loop_edge_count(&multi, loops));
        } else {
            other = igraph_sparsemat_get(adjmatrix, to, from);
            multi = multi > other ? multi : other;
        }
        for (igraph_int_t count = 0; count < multi; count++) {
            IGRAPH_CHECK(igraph_vector_int_push_back(edges, from));
            IGRAPH_CHECK(igraph_vector_int_push_back(edges, to));
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_sparse_adjacency_min(
    igraph_sparsemat_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_loops_t loops
) {
    igraph_sparsemat_iterator_t it;
    igraph_real_t other;

    igraph_sparsemat_iterator_init(&it, adjmatrix);
    for (; !igraph_sparsemat_iterator_end(&it); igraph_sparsemat_iterator_next(&it)) {
        igraph_int_t from = igraph_sparsemat_iterator_row(&it);
        igraph_int_t to = igraph_sparsemat_iterator_col(&it);
        if (to < from) {
            continue;
        }
        igraph_int_t multi = igraph_sparsemat_iterator_get(&it);
        if (to == from) {
            IGRAPH_CHECK(igraph_i_adjust_loop_edge_count(&multi, loops));
        } else {
            other = igraph_sparsemat_get(adjmatrix, to, from);
            multi = multi < other ? multi : other;
        }
        for (igraph_int_t count = 0; count < multi; count++) {
            IGRAPH_CHECK(igraph_vector_int_push_back(edges, from));
            IGRAPH_CHECK(igraph_vector_int_push_back(edges, to));
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_sparse_adjacency_upper(
    igraph_sparsemat_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_loops_t loops
) {
    igraph_sparsemat_iterator_t it;

    /* IGRAPH_LOOPS_TWICE is treated as IGRAPH_LOOPS_ONCE -- it makes no sense
     * for loops to appear twice in the adjacency matrix when the lower triangle
     * is empty; double-counting of loops in undirected graphs happens because
     * the upper and the lower triangle are added on top of each other on the
     * diagonal. See discussion in #2501:
     *
     * https://github.com/igraph/igraph/issues/2501#issuecomment-1949345675 */
    if (loops == IGRAPH_LOOPS_TWICE) {
        loops = IGRAPH_LOOPS_ONCE;
    }

    igraph_sparsemat_iterator_init(&it, adjmatrix);
    for (; !igraph_sparsemat_iterator_end(&it); igraph_sparsemat_iterator_next(&it)) {
        igraph_int_t from = igraph_sparsemat_iterator_row(&it);
        igraph_int_t to = igraph_sparsemat_iterator_col(&it);
        if (to < from) {
            continue;
        }
        igraph_int_t multi = igraph_sparsemat_iterator_get(&it);
        if (to == from) {
            IGRAPH_CHECK(igraph_i_adjust_loop_edge_count(&multi, loops));
        }
        for (igraph_int_t count = 0; count < multi; count++) {
            IGRAPH_CHECK(igraph_vector_int_push_back(edges, from));
            IGRAPH_CHECK(igraph_vector_int_push_back(edges, to));
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_sparse_adjacency_lower(
    igraph_sparsemat_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_loops_t loops
) {
    igraph_sparsemat_iterator_t it;

    /* IGRAPH_LOOPS_TWICE is treated as IGRAPH_LOOPS_ONCE -- it makes no sense
     * for loops to appear twice in the adjacency matrix when the lower triangle
     * is empty; double-counting of loops in undirected graphs happens because
     * the upper and the lower triangle are added on top of each other on the
     * diagonal. See discussion in #2501:
     *
     * https://github.com/igraph/igraph/issues/2501#issuecomment-1949345675 */
    if (loops == IGRAPH_LOOPS_TWICE) {
        loops = IGRAPH_LOOPS_ONCE;
    }

    igraph_sparsemat_iterator_init(&it, adjmatrix);
    for (; !igraph_sparsemat_iterator_end(&it); igraph_sparsemat_iterator_next(&it)) {
        igraph_int_t from = igraph_sparsemat_iterator_row(&it);
        igraph_int_t to = igraph_sparsemat_iterator_col(&it);
        if (to > from) {
            continue;
        }
        igraph_int_t multi = igraph_sparsemat_iterator_get(&it);
        if (to == from) {
            IGRAPH_CHECK(igraph_i_adjust_loop_edge_count(&multi, loops));
        }
        for (igraph_int_t count = 0; count < multi; count++) {
            IGRAPH_CHECK(igraph_vector_int_push_back(edges, from));
            IGRAPH_CHECK(igraph_vector_int_push_back(edges, to));
        }
    }
    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_sparse_adjacency_undirected(
    igraph_sparsemat_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_loops_t loops
) {
    igraph_bool_t sym;

    IGRAPH_CHECK(igraph_sparsemat_is_symmetric(adjmatrix, &sym));
    if (!sym) {
        IGRAPH_ERROR(
            "Adjacency matrix should be symmetric to produce an undirected graph.",
            IGRAPH_EINVAL
        );
    }
    return igraph_i_sparse_adjacency_max(adjmatrix, edges, loops);
}

/**
 * \ingroup generators
 * \function igraph_sparse_adjacency
 * \brief Creates a graph from a sparse adjacency matrix.
 *
 * This has the same functionality as \ref igraph_adjacency(), but uses
 * a column-compressed adjacency matrix.
 *
 * </para><para>
 * Time complexity: O(|E|),
 * where |E| is the number of edges in the graph.
 */

igraph_error_t igraph_sparse_adjacency(igraph_t *graph, igraph_sparsemat_t *adjmatrix,
        igraph_adjacency_t mode, igraph_loops_t loops) {

    igraph_vector_int_t edges = IGRAPH_VECTOR_NULL;
    igraph_int_t no_of_nodes = igraph_sparsemat_nrow(adjmatrix);
    igraph_int_t no_of_nonzeros = igraph_sparsemat_count_nonzero(adjmatrix);
    igraph_int_t approx_no_of_edges;

    if (!igraph_sparsemat_is_cc(adjmatrix)) {
        IGRAPH_ERROR("Sparse adjacency matrix should be in column-compressed "
               "form.", IGRAPH_EINVAL);
    }
    if (no_of_nodes != igraph_sparsemat_ncol(adjmatrix)) {
        IGRAPH_ERROR("Adjacency matrix is non-square.", IGRAPH_EINVAL);
    }

    if (no_of_nodes != 0 && igraph_sparsemat_min(adjmatrix) < 0) {
        IGRAPH_ERRORF("Edge counts should be non-negative, found %g.", IGRAPH_EINVAL,
                igraph_sparsemat_min(adjmatrix));
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

    /* Approximate the number of edges in the graph based on the number of
     * nonzero elements in the matrix */
    switch (mode) {
        case IGRAPH_ADJ_DIRECTED:
        case IGRAPH_ADJ_PLUS:
        case IGRAPH_ADJ_UPPER:
        case IGRAPH_ADJ_LOWER:
            approx_no_of_edges = no_of_nonzeros;
            break;
        case IGRAPH_ADJ_UNDIRECTED:
        case IGRAPH_ADJ_MAX:
        case IGRAPH_ADJ_MIN:
            approx_no_of_edges = no_of_nonzeros / 2;
            break;
        default:
            approx_no_of_edges = no_of_nonzeros;
            break;
    }

    IGRAPH_CHECK(igraph_vector_int_reserve(&edges, approx_no_of_edges * 2));

    /* Collect the edges */
    switch (mode) {
    case IGRAPH_ADJ_DIRECTED:
    case IGRAPH_ADJ_PLUS:
        IGRAPH_CHECK(igraph_i_sparse_adjacency_directed_or_plus(adjmatrix, &edges, mode, loops));
        break;
    case IGRAPH_ADJ_MAX:
        IGRAPH_CHECK(igraph_i_sparse_adjacency_max(adjmatrix, &edges, loops));
        break;
    case IGRAPH_ADJ_UNDIRECTED:
        IGRAPH_CHECK(igraph_i_sparse_adjacency_undirected(adjmatrix, &edges, loops));
        break;
    case IGRAPH_ADJ_UPPER:
        IGRAPH_CHECK(igraph_i_sparse_adjacency_upper(adjmatrix, &edges, loops));
        break;
    case IGRAPH_ADJ_LOWER:
        IGRAPH_CHECK(igraph_i_sparse_adjacency_lower(adjmatrix, &edges, loops));
        break;
    case IGRAPH_ADJ_MIN:
        IGRAPH_CHECK(igraph_i_sparse_adjacency_min(adjmatrix, &edges, loops));
        break;
    default:
        IGRAPH_ERROR("Invalid adjacency mode.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, (mode == IGRAPH_ADJ_DIRECTED)));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_sparse_weighted_adjacency_max (
    igraph_sparsemat_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_vector_t *weights, igraph_loops_t loops
) {
    igraph_sparsemat_iterator_t it;
    igraph_sparsemat_iterator_init(&it, adjmatrix);
    igraph_int_t e = 0;
    igraph_real_t other;

    for (; !igraph_sparsemat_iterator_end(&it); igraph_sparsemat_iterator_next(&it)) {
        igraph_int_t from = igraph_sparsemat_iterator_row(&it);
        igraph_int_t to = igraph_sparsemat_iterator_col(&it);
        if (to < from) {
            continue;
        }
        igraph_real_t weight = igraph_sparsemat_iterator_get(&it);
        if (to == from) {
            igraph_i_adjust_loop_edge_weight(&weight, loops);
        } else {
            other = igraph_sparsemat_get(adjmatrix, to, from);
            weight = weight > other ? weight : other;
        }
        if (weight != 0) {
            VECTOR(*weights)[e/2] = weight;
            VECTOR(*edges)[e++] = from;
            VECTOR(*edges)[e++] = to;
        }
    }
    igraph_vector_int_resize(edges, e); /* shrinks */
    igraph_vector_resize(weights, e/2); /* shrinks */

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_sparse_weighted_adjacency_min (
    igraph_sparsemat_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_vector_t *weights, igraph_loops_t loops
) {
    igraph_sparsemat_iterator_t it;
    igraph_int_t e = 0;
    igraph_real_t other;

    igraph_sparsemat_iterator_init(&it, adjmatrix);
    for (; !igraph_sparsemat_iterator_end(&it); igraph_sparsemat_iterator_next(&it)) {
        igraph_int_t from = igraph_sparsemat_iterator_row(&it);
        igraph_int_t to = igraph_sparsemat_iterator_col(&it);
        if (to < from) {
            continue;
        }
        igraph_real_t weight = igraph_sparsemat_iterator_get(&it);
        if (to == from) {
            igraph_i_adjust_loop_edge_weight(&weight, loops);
        } else {
            other = igraph_sparsemat_get(adjmatrix, to, from);
            weight = weight < other ? weight : other;
        }
        if (weight != 0) {
            VECTOR(*weights)[e/2] = weight;
            VECTOR(*edges)[e++] = from;
            VECTOR(*edges)[e++] = to;
        }
    }
    igraph_vector_int_resize(edges, e); /* shrinks */
    igraph_vector_resize(weights, e/2); /* shrinks */

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_sparse_weighted_adjacency_plus (
    igraph_sparsemat_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_vector_t *weights, igraph_loops_t loops
) {
    igraph_sparsemat_iterator_t it;
    igraph_int_t e = 0;
    igraph_real_t other;

    igraph_sparsemat_iterator_init(&it, adjmatrix);
    for (; !igraph_sparsemat_iterator_end(&it); igraph_sparsemat_iterator_next(&it)) {
        igraph_int_t from = igraph_sparsemat_iterator_row(&it);
        igraph_int_t to = igraph_sparsemat_iterator_col(&it);
        if (to < from) {
            continue;
        }
        igraph_real_t weight = igraph_sparsemat_iterator_get(&it);
        if (to == from) {
            igraph_i_adjust_loop_edge_weight(&weight, loops);
        } else {
            other = igraph_sparsemat_get(adjmatrix, to, from);
            weight += other;
        }
        if (weight != 0) {
            VECTOR(*weights)[e/2] = weight;
            VECTOR(*edges)[e++] = from;
            VECTOR(*edges)[e++] = to;
        }
    }
    igraph_vector_int_resize(edges, e); /* shrinks */
    igraph_vector_resize(weights, e/2); /* shrinks */

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_sparse_weighted_adjacency_upper(
    igraph_sparsemat_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_vector_t *weights, igraph_loops_t loops
) {
    igraph_sparsemat_iterator_t it;
    igraph_sparsemat_iterator_init(&it, adjmatrix);
    igraph_int_t e = 0;

    /* IGRAPH_LOOPS_TWICE is treated as IGRAPH_LOOPS_ONCE -- it makes no sense
     * for loops to appear twice in the adjacency matrix when the lower triangle
     * is empty; double-counting of loops in undirected graphs happens because
     * the upper and the lower triangle are added on top of each other on the
     * diagonal. See discussion in #2501:
     *
     * https://github.com/igraph/igraph/issues/2501#issuecomment-1949345675 */
    if (loops == IGRAPH_LOOPS_TWICE) {
        loops = IGRAPH_LOOPS_ONCE;
    }

    for (; !igraph_sparsemat_iterator_end(&it); igraph_sparsemat_iterator_next(&it)) {
        igraph_int_t from = igraph_sparsemat_iterator_row(&it);
        igraph_int_t to = igraph_sparsemat_iterator_col(&it);
        igraph_real_t weight = igraph_sparsemat_iterator_get(&it);
        if (to < from) {
            continue;
        }
        if (to == from) {
            igraph_i_adjust_loop_edge_weight(&weight, loops);
        }
        if (weight != 0) {
            VECTOR(*weights)[e/2] = weight;
            VECTOR(*edges)[e++] = from;
            VECTOR(*edges)[e++] = to;
        }
    }
    igraph_vector_int_resize(edges, e); /* shrinks */
    igraph_vector_resize(weights, e/2); /* shrinks */

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_sparse_weighted_adjacency_lower(
    igraph_sparsemat_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_vector_t *weights, igraph_loops_t loops
) {
    igraph_sparsemat_iterator_t it;
    igraph_sparsemat_iterator_init(&it, adjmatrix);
    igraph_int_t e = 0;

    /* IGRAPH_LOOPS_TWICE is treated as IGRAPH_LOOPS_ONCE -- it makes no sense
     * for loops to appear twice in the adjacency matrix when the lower triangle
     * is empty; double-counting of loops in undirected graphs happens because
     * the upper and the lower triangle are added on top of each other on the
     * diagonal. See discussion in #2501:
     *
     * https://github.com/igraph/igraph/issues/2501#issuecomment-1949345675 */
    if (loops == IGRAPH_LOOPS_TWICE) {
        loops = IGRAPH_LOOPS_ONCE;
    }

    for (; !igraph_sparsemat_iterator_end(&it); igraph_sparsemat_iterator_next(&it)) {
        igraph_int_t from = igraph_sparsemat_iterator_row(&it);
        igraph_int_t to = igraph_sparsemat_iterator_col(&it);
        igraph_real_t weight = igraph_sparsemat_iterator_get(&it);
        if (to > from) {
            continue;
        }
        if (to == from) {
            igraph_i_adjust_loop_edge_weight(&weight, loops);
        }
        if (weight != 0) {
            VECTOR(*weights)[e/2] = weight;
            VECTOR(*edges)[e++] = from;
            VECTOR(*edges)[e++] = to;
        }
    }
    igraph_vector_int_resize(edges, e); /* shrinks */
    igraph_vector_resize(weights, e/2); /* shrinks */

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_sparse_weighted_adjacency_undirected (
    igraph_sparsemat_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_vector_t *weights, igraph_loops_t loops
) {
    igraph_bool_t sym;

    IGRAPH_CHECK(igraph_sparsemat_is_symmetric(adjmatrix, &sym));
    if (!sym) {
        IGRAPH_ERROR(
            "Adjacency matrix should be symmetric to produce an undirected graph.",
            IGRAPH_EINVAL
        );
    }
    return igraph_i_sparse_weighted_adjacency_max(adjmatrix, edges, weights, loops);
}


static igraph_error_t igraph_i_sparse_weighted_adjacency_directed(
    igraph_sparsemat_t *adjmatrix, igraph_vector_int_t *edges,
    igraph_vector_t *weights, igraph_loops_t loops
) {
    igraph_sparsemat_iterator_t it;
    igraph_sparsemat_iterator_init(&it, adjmatrix);
    igraph_int_t e = 0;

    /* For sake of consistency with the rest of the library, IGRAPH_LOOPS_TWICE
     * is treated as IGRAPH_LOOPS_ONCE for directed graphs */
    if (loops == IGRAPH_LOOPS_TWICE) {
        loops = IGRAPH_LOOPS_ONCE;
    }

    for (; !igraph_sparsemat_iterator_end(&it); igraph_sparsemat_iterator_next(&it)) {
        igraph_int_t from = igraph_sparsemat_iterator_row(&it);
        igraph_int_t to = igraph_sparsemat_iterator_col(&it);
        igraph_real_t weight = igraph_sparsemat_iterator_get(&it);
        if (to == from) {
            igraph_i_adjust_loop_edge_weight(&weight, loops);
        }
        if (weight != 0) {
            VECTOR(*weights)[e/2] = weight;
            VECTOR(*edges)[e++] = from;
            VECTOR(*edges)[e++] = to;
        }
    }
    igraph_vector_int_resize(edges, e); /* shrinks */
    igraph_vector_resize(weights, e/2); /* shrinks */

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup generators
 * \function igraph_sparse_weighted_adjacency
 * \brief Creates a graph from a weighted sparse adjacency matrix.
 *
 * This has the same functionality as \ref igraph_weighted_adjacency(), but uses
 * a column-compressed adjacency matrix.
 *
 * </para><para>
 * Time complexity: O(|E|),
 * where |E| is the number of edges in the graph.
 */


igraph_error_t igraph_sparse_weighted_adjacency(
    igraph_t *graph, igraph_sparsemat_t *adjmatrix, igraph_adjacency_t mode,
    igraph_vector_t *weights, igraph_loops_t loops
) {
    igraph_vector_int_t edges;
    igraph_int_t no_of_nodes = igraph_sparsemat_nrow(adjmatrix);
    igraph_int_t no_of_edges = igraph_sparsemat_count_nonzero(adjmatrix);

    if (!igraph_sparsemat_is_cc(adjmatrix)) {
        IGRAPH_ERROR("Sparse adjacency matrix should be in column-compressed form.", IGRAPH_EINVAL);
    }
    if (no_of_nodes != igraph_sparsemat_ncol(adjmatrix)) {
        IGRAPH_ERROR("Adjacency matrix is non-square.", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges * 2);
    IGRAPH_CHECK(igraph_vector_resize(weights, no_of_edges));

    /* Collect the edges */
    switch (mode) {
    case IGRAPH_ADJ_DIRECTED:
        IGRAPH_CHECK(igraph_i_sparse_weighted_adjacency_directed(adjmatrix, &edges,
                     weights, loops));
        break;
    case IGRAPH_ADJ_MAX:
        IGRAPH_CHECK(igraph_i_sparse_weighted_adjacency_max(adjmatrix, &edges,
                     weights, loops));
        break;
    case IGRAPH_ADJ_UNDIRECTED:
        IGRAPH_CHECK(igraph_i_sparse_weighted_adjacency_undirected(adjmatrix, &edges,
                     weights, loops));
        break;
    case IGRAPH_ADJ_UPPER:
        IGRAPH_CHECK(igraph_i_sparse_weighted_adjacency_upper(adjmatrix, &edges,
                     weights, loops));
        break;
    case IGRAPH_ADJ_LOWER:
        IGRAPH_CHECK(igraph_i_sparse_weighted_adjacency_lower(adjmatrix, &edges,
                     weights, loops));
        break;
    case IGRAPH_ADJ_MIN:
        IGRAPH_CHECK(igraph_i_sparse_weighted_adjacency_min(adjmatrix, &edges,
                     weights, loops));
        break;
    case IGRAPH_ADJ_PLUS:
        IGRAPH_CHECK(igraph_i_sparse_weighted_adjacency_plus(adjmatrix, &edges,
                     weights, loops));
        break;
    default:
        IGRAPH_ERROR("Invalid adjacency mode.", IGRAPH_EINVAL);
    }

    /* Create graph */
    IGRAPH_CHECK(igraph_empty(graph, no_of_nodes, (mode == IGRAPH_ADJ_DIRECTED)));
    IGRAPH_FINALLY(igraph_destroy, graph);
    if (igraph_vector_int_size(&edges) > 0) {
        IGRAPH_CHECK(igraph_add_edges(graph, &edges, NULL));
    }
    IGRAPH_FINALLY_CLEAN(1);

    /* Cleanup */
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
