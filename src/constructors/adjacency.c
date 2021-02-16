/* -*- mode: C -*-  */
/*
   IGraph library.
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
#include "igraph_attributes.h"
#include "igraph_interface.h"

static int igraph_i_adjacency_directed(igraph_matrix_t *adjmatrix,
                                       igraph_vector_t *edges);
static int igraph_i_adjacency_max(igraph_matrix_t *adjmatrix,
                                  igraph_vector_t *edges);
static int igraph_i_adjacency_upper(igraph_matrix_t *adjmatrix,
                                    igraph_vector_t *edges);
static int igraph_i_adjacency_lower(igraph_matrix_t *adjmatrix,
                                    igraph_vector_t *edges);
static int igraph_i_adjacency_min(igraph_matrix_t *adjmatrix,
                                  igraph_vector_t *edges);

static int igraph_i_adjacency_directed(igraph_matrix_t *adjmatrix, igraph_vector_t *edges) {

    long int no_of_nodes = igraph_matrix_nrow(adjmatrix);
    long int i, j, k;

    for (i = 0; i < no_of_nodes; i++) {
        for (j = 0; j < no_of_nodes; j++) {
            long int M = (long int) MATRIX(*adjmatrix, i, j);
            for (k = 0; k < M; k++) {
                IGRAPH_CHECK(igraph_vector_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(edges, j));
            }
        }
    }

    return 0;
}

static int igraph_i_adjacency_max(igraph_matrix_t *adjmatrix, igraph_vector_t *edges) {

    long int no_of_nodes = igraph_matrix_nrow(adjmatrix);
    long int i, j, k;

    for (i = 0; i < no_of_nodes; i++) {
        for (j = i; j < no_of_nodes; j++) {
            long int M1 = (long int) MATRIX(*adjmatrix, i, j);
            long int M2 = (long int) MATRIX(*adjmatrix, j, i);
            if (M1 < M2) {
                M1 = M2;
            }
            for (k = 0; k < M1; k++) {
                IGRAPH_CHECK(igraph_vector_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(edges, j));
            }
        }
    }

    return 0;
}

static int igraph_i_adjacency_upper(igraph_matrix_t *adjmatrix, igraph_vector_t *edges) {

    long int no_of_nodes = igraph_matrix_nrow(adjmatrix);
    long int i, j, k;

    for (i = 0; i < no_of_nodes; i++) {
        for (j = i; j < no_of_nodes; j++) {
            long int M = (long int) MATRIX(*adjmatrix, i, j);
            for (k = 0; k < M; k++) {
                IGRAPH_CHECK(igraph_vector_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(edges, j));
            }
        }
    }
    return 0;
}

static int igraph_i_adjacency_lower(igraph_matrix_t *adjmatrix, igraph_vector_t *edges) {

    long int no_of_nodes = igraph_matrix_nrow(adjmatrix);
    long int i, j, k;

    for (i = 0; i < no_of_nodes; i++) {
        for (j = 0; j <= i; j++) {
            long int M = (long int) MATRIX(*adjmatrix, i, j);
            for (k = 0; k < M; k++) {
                IGRAPH_CHECK(igraph_vector_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(edges, j));
            }
        }
    }
    return 0;
}

static int igraph_i_adjacency_min(igraph_matrix_t *adjmatrix, igraph_vector_t *edges) {

    long int no_of_nodes = igraph_matrix_nrow(adjmatrix);
    long int i, j, k;

    for (i = 0; i < no_of_nodes; i++) {
        for (j = i; j < no_of_nodes; j++) {
            long int M1 = (long int) MATRIX(*adjmatrix, i, j);
            long int M2 = (long int) MATRIX(*adjmatrix, j, i);
            if (M1 > M2) {
                M1 = M2;
            }
            for (k = 0; k < M1; k++) {
                IGRAPH_CHECK(igraph_vector_push_back(edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(edges, j));
            }
        }
    }

    return 0;
}

/**
 * \ingroup generators
 * \function igraph_adjacency
 * \brief Creates a graph from an adjacency matrix.
 *
 * The order of the vertices in the matrix is preserved, i.e. the vertex
 * corresponding to the first row/column will be vertex with id 0, the
 * next row is for vertex 1, etc.
 * \param graph Pointer to an uninitialized graph object.
 * \param adjmatrix The adjacency matrix. How it is interpreted
 *        depends on the \p mode argument.
 * \param mode Constant to specify how the given matrix is interpreted
 *        as an adjacency matrix. Possible values
 *        (A(i,j)
 *        is the element in row i and column
 *        j in the adjacency matrix
 *        \p adjmatrix):
 *        \clist
 *        \cli IGRAPH_ADJ_DIRECTED
 *          the graph will be directed and
 *          an element gives the number of edges between two vertices.
 *        \cli IGRAPH_ADJ_UNDIRECTED
 *          this is the same as \c IGRAPH_ADJ_MAX,
 *          for convenience.
 *        \cli IGRAPH_ADJ_MAX
 *          undirected graph will be created
 *          and the number of edges between vertices
 *          i and
 *          j is
 *          max(A(i,j), A(j,i)).
 *        \cli IGRAPH_ADJ_MIN
 *          undirected graph will be created
 *          with min(A(i,j), A(j,i))
 *          edges between vertices
 *          i and
 *          j.
 *        \cli IGRAPH_ADJ_PLUS
 *          undirected graph will be created
 *          with A(i,j)+A(j,i) edges
 *          between vertices
 *          i and
 *          j.
 *        \cli IGRAPH_ADJ_UPPER
 *          undirected graph will be created,
 *          only the upper right triangle (including the diagonal) is
 *          used for the number of edges.
 *        \cli IGRAPH_ADJ_LOWER
 *          undirected graph will be created,
 *          only the lower left triangle (including the diagonal) is
 *          used for creating the edges.
 *       \endclist
 * \return Error code,
 *         \c IGRAPH_NONSQUARE: non-square matrix.
 *
 * Time complexity: O(|V||V|),
 * |V| is the number of vertices in the graph.
 *
 * \example examples/simple/igraph_adjacency.c
 */
int igraph_adjacency(igraph_t *graph, igraph_matrix_t *adjmatrix,
                     igraph_adjacency_t mode) {

    igraph_vector_t edges = IGRAPH_VECTOR_NULL;
    long int no_of_nodes;

    /* Some checks */
    if (igraph_matrix_nrow(adjmatrix) != igraph_matrix_ncol(adjmatrix)) {
        IGRAPH_ERROR("Non-square matrix", IGRAPH_NONSQUARE);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

    /* Collect the edges */
    no_of_nodes = igraph_matrix_nrow(adjmatrix);
    switch (mode) {
    case IGRAPH_ADJ_DIRECTED:
        IGRAPH_CHECK(igraph_i_adjacency_directed(adjmatrix, &edges));
        break;
    case IGRAPH_ADJ_MAX:
        IGRAPH_CHECK(igraph_i_adjacency_max(adjmatrix, &edges));
        break;
    case IGRAPH_ADJ_UPPER:
        IGRAPH_CHECK(igraph_i_adjacency_upper(adjmatrix, &edges));
        break;
    case IGRAPH_ADJ_LOWER:
        IGRAPH_CHECK(igraph_i_adjacency_lower(adjmatrix, &edges));
        break;
    case IGRAPH_ADJ_MIN:
        IGRAPH_CHECK(igraph_i_adjacency_min(adjmatrix, &edges));
        break;
    case IGRAPH_ADJ_PLUS:
        IGRAPH_CHECK(igraph_i_adjacency_directed(adjmatrix, &edges));
        break;
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, (igraph_integer_t) no_of_nodes,
                               (mode == IGRAPH_ADJ_DIRECTED)));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

static int igraph_i_weighted_adjacency_directed(
        const igraph_matrix_t *adjmatrix,
        igraph_vector_t *edges,
        igraph_vector_t *weights,
        igraph_bool_t loops);
static int igraph_i_weighted_adjacency_plus(
        const igraph_matrix_t *adjmatrix,
        igraph_vector_t *edges,
        igraph_vector_t *weights,
        igraph_bool_t loops);
static int igraph_i_weighted_adjacency_max(
        const igraph_matrix_t *adjmatrix,
        igraph_vector_t *edges,
        igraph_vector_t *weights,
        igraph_bool_t loops);
static int igraph_i_weighted_adjacency_upper(
        const igraph_matrix_t *adjmatrix,
        igraph_vector_t *edges,
        igraph_vector_t *weights,
        igraph_bool_t loops);
static int igraph_i_weighted_adjacency_lower(
        const igraph_matrix_t *adjmatrix,
        igraph_vector_t *edges,
        igraph_vector_t *weights,
        igraph_bool_t loops);
static int igraph_i_weighted_adjacency_min(
        const igraph_matrix_t *adjmatrix,
        igraph_vector_t *edges,
        igraph_vector_t *weights,
        igraph_bool_t loops);

static int igraph_i_weighted_adjacency_directed(
        const igraph_matrix_t *adjmatrix,
        igraph_vector_t *edges,
        igraph_vector_t *weights,
        igraph_bool_t loops) {

    long int no_of_nodes = igraph_matrix_nrow(adjmatrix);
    long int i, j;

    for (i = 0; i < no_of_nodes; i++) {
        for (j = 0; j < no_of_nodes; j++) {
            igraph_real_t M = MATRIX(*adjmatrix, i, j);
            if (M == 0.0) {
                continue;
            }
            if (i == j && !loops) {
                continue;
            }
            IGRAPH_CHECK(igraph_vector_push_back(edges, i));
            IGRAPH_CHECK(igraph_vector_push_back(edges, j));
            IGRAPH_CHECK(igraph_vector_push_back(weights, M));
        }
    }

    return 0;
}

static int igraph_i_weighted_adjacency_plus(
        const igraph_matrix_t *adjmatrix,
        igraph_vector_t *edges,
        igraph_vector_t *weights,
        igraph_bool_t loops) {

    long int no_of_nodes = igraph_matrix_nrow(adjmatrix);
    long int i, j;

    for (i = 0; i < no_of_nodes; i++) {
        for (j = i; j < no_of_nodes; j++) {
            igraph_real_t M = MATRIX(*adjmatrix, i, j) + MATRIX(*adjmatrix, j, i);
            if (M == 0.0) {
                continue;
            }
            if (i == j && !loops) {
                continue;
            }
            if (i == j) {
                M /= 2;
            }
            IGRAPH_CHECK(igraph_vector_push_back(edges, i));
            IGRAPH_CHECK(igraph_vector_push_back(edges, j));
            IGRAPH_CHECK(igraph_vector_push_back(weights, M));
        }
    }

    return 0;
}

static int igraph_i_weighted_adjacency_max(
        const igraph_matrix_t *adjmatrix,
        igraph_vector_t *edges,
        igraph_vector_t *weights,
        igraph_bool_t loops) {

    long int no_of_nodes = igraph_matrix_nrow(adjmatrix);
    long int i, j;

    for (i = 0; i < no_of_nodes; i++) {
        for (j = i; j < no_of_nodes; j++) {
            igraph_real_t M1 = MATRIX(*adjmatrix, i, j);
            igraph_real_t M2 = MATRIX(*adjmatrix, j, i);
            if (M1 < M2) {
                M1 = M2;
            }
            if (M1 == 0.0) {
                continue;
            }
            if (i == j && !loops) {
                continue;
            }
            IGRAPH_CHECK(igraph_vector_push_back(edges, i));
            IGRAPH_CHECK(igraph_vector_push_back(edges, j));
            IGRAPH_CHECK(igraph_vector_push_back(weights, M1));
        }
    }
    return 0;
}

static int igraph_i_weighted_adjacency_upper(
        const igraph_matrix_t *adjmatrix,
        igraph_vector_t *edges,
        igraph_vector_t *weights,
        igraph_bool_t loops) {

    long int no_of_nodes = igraph_matrix_nrow(adjmatrix);
    long int i, j;

    for (i = 0; i < no_of_nodes; i++) {
        for (j = i; j < no_of_nodes; j++) {
            igraph_real_t M = MATRIX(*adjmatrix, i, j);
            if (M == 0.0) {
                continue;
            }
            if (i == j && !loops) {
                continue;
            }
            IGRAPH_CHECK(igraph_vector_push_back(edges, i));
            IGRAPH_CHECK(igraph_vector_push_back(edges, j));
            IGRAPH_CHECK(igraph_vector_push_back(weights, M));
        }
    }
    return 0;
}

static int igraph_i_weighted_adjacency_lower(
        const igraph_matrix_t *adjmatrix,
        igraph_vector_t *edges,
        igraph_vector_t *weights,
        igraph_bool_t loops) {

    long int no_of_nodes = igraph_matrix_nrow(adjmatrix);
    long int i, j;

    for (i = 0; i < no_of_nodes; i++) {
        for (j = 0; j <= i; j++) {
            igraph_real_t M = MATRIX(*adjmatrix, i, j);
            if (M == 0.0) {
                continue;
            }
            if (i == j && !loops) {
                continue;
            }
            IGRAPH_CHECK(igraph_vector_push_back(edges, i));
            IGRAPH_CHECK(igraph_vector_push_back(edges, j));
            IGRAPH_CHECK(igraph_vector_push_back(weights, M));
        }
    }
    return 0;
}

static int igraph_i_weighted_adjacency_min(
        const igraph_matrix_t *adjmatrix,
        igraph_vector_t *edges,
        igraph_vector_t *weights,
        igraph_bool_t loops) {

    long int no_of_nodes = igraph_matrix_nrow(adjmatrix);
    long int i, j;

    for (i = 0; i < no_of_nodes; i++) {
        for (j = i; j < no_of_nodes; j++) {
            igraph_real_t M1 = MATRIX(*adjmatrix, i, j);
            igraph_real_t M2 = MATRIX(*adjmatrix, j, i);
            if (M1 > M2) {
                M1 = M2;
            }
            if (M1 == 0.0) {
                continue;
            }
            if (i == j && !loops) {
                continue;
            }
            IGRAPH_CHECK(igraph_vector_push_back(edges, i));
            IGRAPH_CHECK(igraph_vector_push_back(edges, j));
            IGRAPH_CHECK(igraph_vector_push_back(weights, M1));
        }
    }

    return 0;
}

/**
 * \ingroup generators
 * \function igraph_weighted_adjacency
 * \brief Creates a graph from a weighted adjacency matrix.
 *
 * The order of the vertices in the matrix is preserved, i.e. the vertex
 * corresponding to the first row/column will be vertex with id 0, the
 * next row is for vertex 1, etc.
 * \param graph Pointer to an uninitialized graph object.
 * \param adjmatrix The weighted adjacency matrix. How it is interpreted
 *        depends on the \p mode argument. The common feature is that
 *        edges with zero weights are considered nonexistent (however,
 *        negative weights are permitted).
 * \param mode Constant to specify how the given matrix is interpreted
 *        as an adjacency matrix. Possible values
 *        (A(i,j)
 *        is the element in row i and column
 *        j in the adjacency matrix
 *        \p adjmatrix):
 *        \clist
 *        \cli IGRAPH_ADJ_DIRECTED
 *          the graph will be directed and
 *          an element gives the weight of the edge between two vertices.
 *        \cli IGRAPH_ADJ_UNDIRECTED
 *          this is the same as \c IGRAPH_ADJ_MAX,
 *          for convenience.
 *        \cli IGRAPH_ADJ_MAX
 *          undirected graph will be created
 *          and the weight of the edge between vertices
 *          i and
 *          j is
 *          max(A(i,j), A(j,i)).
 *        \cli IGRAPH_ADJ_MIN
 *          undirected graph will be created
 *          with edge weight min(A(i,j), A(j,i))
 *          between vertices
 *          i and
 *          j.
 *        \cli IGRAPH_ADJ_PLUS
 *          undirected graph will be created
 *          with edge weight A(i,j)+A(j,i)
 *          between vertices
 *          i and
 *          j.
 *        \cli IGRAPH_ADJ_UPPER
 *          undirected graph will be created,
 *          only the upper right triangle (including the diagonal) is
 *          used for the edge weights.
 *        \cli IGRAPH_ADJ_LOWER
 *          undirected graph will be created,
 *          only the lower left triangle (including the diagonal) is
 *          used for the edge weights.
 *       \endclist
 * \param attr the name of the attribute that will store the edge weights.
 *         If \c NULL , it will use \c weight as the attribute name.
 * \param loops Logical scalar, whether to ignore the diagonal elements
 *         in the adjacency matrix.
 * \return Error code,
 *         \c IGRAPH_NONSQUARE: non-square matrix.
 *
 * Time complexity: O(|V||V|),
 * |V| is the number of vertices in the graph.
 *
 * \example examples/simple/igraph_weighted_adjacency.c
 */
int igraph_weighted_adjacency(igraph_t *graph, igraph_matrix_t *adjmatrix,
                              igraph_adjacency_t mode, const char* attr,
                              igraph_bool_t loops) {

    igraph_vector_t edges = IGRAPH_VECTOR_NULL;
    igraph_vector_t weights = IGRAPH_VECTOR_NULL;
    const char* default_attr = "weight";
    igraph_vector_ptr_t attr_vec;
    igraph_attribute_record_t attr_rec;
    long int no_of_nodes;

    /* Some checks */
    if (igraph_matrix_nrow(adjmatrix) != igraph_matrix_ncol(adjmatrix)) {
        IGRAPH_ERROR("Non-square matrix", IGRAPH_NONSQUARE);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&weights, 0);
    IGRAPH_VECTOR_PTR_INIT_FINALLY(&attr_vec, 1);

    /* Collect the edges */
    no_of_nodes = igraph_matrix_nrow(adjmatrix);
    switch (mode) {
    case IGRAPH_ADJ_DIRECTED:
        IGRAPH_CHECK(igraph_i_weighted_adjacency_directed(adjmatrix, &edges,
                     &weights, loops));
        break;
    case IGRAPH_ADJ_MAX:
        IGRAPH_CHECK(igraph_i_weighted_adjacency_max(adjmatrix, &edges,
                     &weights, loops));
        break;
    case IGRAPH_ADJ_UPPER:
        IGRAPH_CHECK(igraph_i_weighted_adjacency_upper(adjmatrix, &edges,
                     &weights, loops));
        break;
    case IGRAPH_ADJ_LOWER:
        IGRAPH_CHECK(igraph_i_weighted_adjacency_lower(adjmatrix, &edges,
                     &weights, loops));
        break;
    case IGRAPH_ADJ_MIN:
        IGRAPH_CHECK(igraph_i_weighted_adjacency_min(adjmatrix, &edges,
                     &weights, loops));
        break;
    case IGRAPH_ADJ_PLUS:
        IGRAPH_CHECK(igraph_i_weighted_adjacency_plus(adjmatrix, &edges,
                     &weights, loops));
        break;
    }

    /* Prepare attribute record */
    attr_rec.name = attr ? attr : default_attr;
    attr_rec.type = IGRAPH_ATTRIBUTE_NUMERIC;
    attr_rec.value = &weights;
    VECTOR(attr_vec)[0] = &attr_rec;

    /* Create graph */
    IGRAPH_CHECK(igraph_empty(graph, (igraph_integer_t) no_of_nodes,
                              (mode == IGRAPH_ADJ_DIRECTED)));
    IGRAPH_FINALLY(igraph_destroy, graph);
    if (igraph_vector_size(&edges) > 0) {
        IGRAPH_CHECK(igraph_add_edges(graph, &edges, &attr_vec));
    }
    IGRAPH_FINALLY_CLEAN(1);

    /* Cleanup */
    igraph_vector_destroy(&edges);
    igraph_vector_destroy(&weights);
    igraph_vector_ptr_destroy(&attr_vec);
    IGRAPH_FINALLY_CLEAN(3);

    return 0;
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
 * \param duplicate Logical, for undirected graphs this specified
 *        whether each edge is included twice, in the vectors of
 *        both adjacent vertices. If this is false (0), then it is
 *        assumed that every edge is included only once. This argument
 *        is ignored for directed graphs.
 * \return Error code.
 *
 * \sa \ref igraph_adjlist_init() for the opposite operation.
 *
 * Time complexity: O(|V|+|E|).
 *
 */
int igraph_adjlist(igraph_t *graph, const igraph_adjlist_t *adjlist,
                   igraph_neimode_t mode, igraph_bool_t duplicate) {

    long int no_of_nodes = igraph_adjlist_size(adjlist);
    long int no_of_edges = 0;
    long int i;

    igraph_vector_t edges;
    long int edgeptr = 0;

    duplicate = duplicate && (mode == IGRAPH_ALL); /* only duplicate if undirected */

    for (i = 0; i < no_of_nodes; i++) {
        no_of_edges += igraph_vector_int_size(igraph_adjlist_get(adjlist, i));
    }

    if (duplicate) {
        no_of_edges /= 2;
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 2 * no_of_edges);

    for (i = 0; i < no_of_nodes; i++) {
        igraph_vector_int_t *neis = igraph_adjlist_get(adjlist, i);
        long int j, n = igraph_vector_int_size(neis);
        long int loops = 0;

        for (j = 0; j < n; j++) {
            long int nei = (long int) VECTOR(*neis)[j];
            if (nei == i) {
                loops++;
            } else {
                if (! duplicate || nei > i) {
                    if (edgeptr + 2 > 2 * no_of_edges) {
                        IGRAPH_ERROR("Invalid adjacency list, most probably not correctly"
                                     " duplicated edges for an undirected graph", IGRAPH_EINVAL);
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
                         " duplicated edges for an undirected graph", IGRAPH_EINVAL);
        }
        for (j = 0; j < loops; j++) {
            VECTOR(edges)[edgeptr++] = i;
            VECTOR(edges)[edgeptr++] = i;
        }
    }

    if (mode == IGRAPH_ALL)
        IGRAPH_CHECK(igraph_create(graph, &edges,
                                   (igraph_integer_t) no_of_nodes, 0));
    else
        IGRAPH_CHECK(igraph_create(graph, &edges,
                                   (igraph_integer_t) no_of_nodes, 1));

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}
