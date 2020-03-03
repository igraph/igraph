/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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
#include "igraph_structural.h"
#include "igraph_memory.h"
#include "igraph_interface.h"
#include "igraph_attributes.h"
#include "igraph_adjlist.h"
#include "igraph_interrupt_internal.h"
#include "igraph_dqueue.h"
#include "config.h"

#include <stdarg.h>
#include <math.h>
#include <string.h>

/**
 * \section about_generators
 *
 * <para>Graph generators create graphs.</para>
 *
 * <para>Almost all functions which create graph objects are documented
 * here. The exceptions are \ref igraph_subgraph() and alike, these
 * create graphs based on another graph.</para>
 */


/**
 * \ingroup generators
 * \function igraph_create
 * \brief Creates a graph with the specified edges.
 *
 * \param graph An uninitialized graph object.
 * \param edges The edges to add, the first two elements are the first
 *        edge, etc.
 * \param n The number of vertices in the graph, if smaller or equal
 *        to the highest vertex id in the \p edges vector it
 *        will be increased automatically. So it is safe to give 0
 *        here.
 * \param directed Boolean, whether to create a directed graph or
 *        not. If yes, then the first edge points from the first
 *        vertex id in \p edges to the second, etc.
 * \return Error code:
 *         \c IGRAPH_EINVEVECTOR: invalid edges
 *         vector (odd number of vertices).
 *         \c IGRAPH_EINVVID: invalid (negative)
 *         vertex id.
 *
 * Time complexity: O(|V|+|E|),
 * |V| is the number of vertices,
 * |E| the number of edges in the
 * graph.
 *
 * \example examples/simple/igraph_create.c
 */
int igraph_create(igraph_t *graph, const igraph_vector_t *edges,
                  igraph_integer_t n, igraph_bool_t directed) {
    igraph_bool_t has_edges = igraph_vector_size(edges) > 0;
    igraph_real_t max = has_edges ? igraph_vector_max(edges) + 1 : 0;

    if (igraph_vector_size(edges) % 2 != 0) {
        IGRAPH_ERROR("Invalid (odd) edges vector", IGRAPH_EINVEVECTOR);
    }
    if (has_edges && !igraph_vector_isininterval(edges, 0, max - 1)) {
        IGRAPH_ERROR("Invalid (negative) vertex id", IGRAPH_EINVVID);
    }

    IGRAPH_CHECK(igraph_empty(graph, n, directed));
    IGRAPH_FINALLY(igraph_destroy, graph);
    if (has_edges) {
        igraph_integer_t vc = igraph_vcount(graph);
        if (vc < max) {
            IGRAPH_CHECK(igraph_add_vertices(graph, (igraph_integer_t) (max - vc), 0));
        }
        IGRAPH_CHECK(igraph_add_edges(graph, edges, 0));
    }

    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

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
 * \brief Creates a graph object from an adjacency matrix.
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
 * \brief Creates a graph object from a weighted adjacency matrix.
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
 * \ingroup generators
 * \function igraph_star
 * \brief Creates a \em star graph, every vertex connects only to the center.
 *
 * \param graph Pointer to an uninitialized graph object, this will
 *        be the result.
 * \param n Integer constant, the number of vertices in the graph.
 * \param mode Constant, gives the type of the star graph to
 *        create. Possible values:
 *        \clist
 *        \cli IGRAPH_STAR_OUT
 *          directed star graph, edges point
 *          \em from the center to the other vertices.
 *        \cli IGRAPH_STAR_IN
 *          directed star graph, edges point
 *          \em to the center from the other vertices.
 *        \cli IGRAPH_STAR_MUTUAL
 *          directed star graph with mutual edges.
 *        \cli IGRAPH_STAR_UNDIRECTED
 *          an undirected star graph is
 *          created.
 *        \endclist
 * \param center Id of the vertex which will be the center of the
 *          graph.
 * \return Error code:
 *         \clist
 *         \cli IGRAPH_EINVVID
 *           invalid number of vertices.
 *         \cli IGRAPH_EINVAL
 *           invalid center vertex.
 *         \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *         \endclist
 *
 * Time complexity: O(|V|), the
 * number of vertices in the graph.
 *
 * \sa \ref igraph_lattice(), \ref igraph_ring(), \ref igraph_tree()
 * for creating other regular structures.
 *
 * \example examples/simple/igraph_star.c
 */

int igraph_star(igraph_t *graph, igraph_integer_t n, igraph_star_mode_t mode,
                igraph_integer_t center) {

    igraph_vector_t edges = IGRAPH_VECTOR_NULL;
    long int i;

    if (n < 0) {
        IGRAPH_ERROR("Invalid number of vertices", IGRAPH_EINVVID);
    }
    if (center < 0 || center > n - 1) {
        IGRAPH_ERROR("Invalid center vertex", IGRAPH_EINVAL);
    }
    if (mode != IGRAPH_STAR_OUT && mode != IGRAPH_STAR_IN &&
        mode != IGRAPH_STAR_MUTUAL && mode != IGRAPH_STAR_UNDIRECTED) {
        IGRAPH_ERROR("invalid mode", IGRAPH_EINVMODE);
    }

    if (mode != IGRAPH_STAR_MUTUAL) {
        IGRAPH_VECTOR_INIT_FINALLY(&edges, (n - 1) * 2);
    } else {
        IGRAPH_VECTOR_INIT_FINALLY(&edges, (n - 1) * 2 * 2);
    }

    if (mode == IGRAPH_STAR_OUT) {
        for (i = 0; i < center; i++) {
            VECTOR(edges)[2 * i] = center;
            VECTOR(edges)[2 * i + 1] = i;
        }
        for (i = center + 1; i < n; i++) {
            VECTOR(edges)[2 * (i - 1)] = center;
            VECTOR(edges)[2 * (i - 1) + 1] = i;
        }
    } else if (mode == IGRAPH_STAR_MUTUAL) {
        for (i = 0; i < center; i++) {
            VECTOR(edges)[4 * i] = center;
            VECTOR(edges)[4 * i + 1] = i;
            VECTOR(edges)[4 * i + 2] = i;
            VECTOR(edges)[4 * i + 3] = center;
        }
        for (i = center + 1; i < n; i++) {
            VECTOR(edges)[4 * i - 4] = center;
            VECTOR(edges)[4 * i - 3] = i;
            VECTOR(edges)[4 * i - 2] = i;
            VECTOR(edges)[4 * i - 1] = center;
        }
    } else {
        for (i = 0; i < center; i++) {
            VECTOR(edges)[2 * i + 1] = center;
            VECTOR(edges)[2 * i] = i;
        }
        for (i = center + 1; i < n; i++) {
            VECTOR(edges)[2 * (i - 1) + 1] = center;
            VECTOR(edges)[2 * (i - 1)] = i;
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, 0,
                               (mode != IGRAPH_STAR_UNDIRECTED)));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \ingroup generators
 * \function igraph_lattice
 * \brief Creates most kinds of lattices.
 *
 * \param graph An uninitialized graph object.
 * \param dimvector Vector giving the sizes of the lattice in each of
 *        its dimensions. Ie. the dimension of the lattice will be the
 *        same as the length of this vector.
 * \param nei Integer value giving the distance (number of steps)
 *        within which two vertices will be connected.
 * \param directed Boolean, whether to create a directed graph. The
 *        direction of the edges is determined by the generation
 *        algorithm and is unlikely to suit you, so this isn't a very
 *        useful option.
 * \param mutual Boolean, if the graph is directed this gives whether
 *        to create all connections as mutual.
 * \param circular Boolean, defines whether the generated lattice is
 *        periodic.
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid (negative)
 *         dimension vector.
 *
 * Time complexity: if \p nei is less than two then it is O(|V|+|E|) (as
 * far as I remember), |V| and |E| are the number of vertices
 * and edges in the generated graph. Otherwise it is O(|V|*d^o+|E|), d
 * is the average degree of the graph, o is the \p nei argument.
 */
int igraph_lattice(igraph_t *graph, const igraph_vector_t *dimvector,
                   igraph_integer_t nei, igraph_bool_t directed, igraph_bool_t mutual,
                   igraph_bool_t circular) {

    long int dims = igraph_vector_size(dimvector);
    long int no_of_nodes = (long int) igraph_vector_prod(dimvector);
    igraph_vector_t edges = IGRAPH_VECTOR_NULL;
    long int *coords, *weights;
    long int i, j;
    int carry, pos;

    if (igraph_vector_any_smaller(dimvector, 0)) {
        IGRAPH_ERROR("Invalid dimension vector", IGRAPH_EINVAL);
    }

    /* init coords & weights */

    coords = igraph_Calloc(dims, long int);
    if (coords == 0) {
        IGRAPH_ERROR("lattice failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, coords);
    weights = igraph_Calloc(dims, long int);
    if (weights == 0) {
        IGRAPH_ERROR("lattice failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, weights);
    if (dims > 0) {
        weights[0] = 1;
        for (i = 1; i < dims; i++) {
            weights[i] = weights[i - 1] * (long int) VECTOR(*dimvector)[i - 1];
        }
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&edges, no_of_nodes * dims +
                                       mutual * directed * no_of_nodes * dims));

    for (i = 0; i < no_of_nodes; i++) {
        IGRAPH_ALLOW_INTERRUPTION();
        for (j = 0; j < dims; j++) {
            if (circular || coords[j] != VECTOR(*dimvector)[j] - 1) {
                long int new_nei;
                if (coords[j] != VECTOR(*dimvector)[j] - 1) {
                    new_nei = i + weights[j] + 1;
                } else {
                    new_nei = i - (long int) (VECTOR(*dimvector)[j] - 1) * weights[j] + 1;
                }
                if (new_nei != i + 1 &&
                    (VECTOR(*dimvector)[j] != 2 || coords[j] != 1 || directed)) {
                    igraph_vector_push_back(&edges, i); /* reserved */
                    igraph_vector_push_back(&edges, new_nei - 1); /* reserved */
                }
            } /* if circular || coords[j] */
            if (mutual && directed && (circular || coords[j] != 0)) {
                long int new_nei;
                if (coords[j] != 0) {
                    new_nei = i - weights[j] + 1;
                } else {
                    new_nei = i + (long int) (VECTOR(*dimvector)[j] - 1) * weights[j] + 1;
                }
                if (new_nei != i + 1 &&
                    (VECTOR(*dimvector)[j] != 2 || !circular)) {
                    igraph_vector_push_back(&edges, i); /* reserved */
                    igraph_vector_push_back(&edges, new_nei - 1); /* reserved */
                }
            } /* if circular || coords[0] */
        } /* for j<dims */

        /* increase coords */
        carry = 1;
        pos = 0;

        while (carry == 1 && pos != dims) {
            if (coords[pos] != VECTOR(*dimvector)[pos] - 1) {
                coords[pos]++;
                carry = 0;
            } else {
                coords[pos] = 0;
                pos++;
            }
        }

    } /* for i<no_of_nodes */

    IGRAPH_CHECK(igraph_create(graph, &edges, (igraph_integer_t) no_of_nodes,
                               directed));
    if (nei >= 2) {
        IGRAPH_CHECK(igraph_connect_neighborhood(graph, nei, IGRAPH_ALL));
    }

    /* clean up */
    igraph_Free(coords);
    igraph_Free(weights);
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(3);

    return 0;
}

/**
 * \ingroup generators
 * \function igraph_ring
 * \brief Creates a \em ring graph, a one dimensional lattice.
 *
 * An undirected (circular) ring on n vertices is commonly known in graph
 * theory as the cycle graph C_n.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param n The number of vertices in the ring.
 * \param directed Logical, whether to create a directed ring.
 * \param mutual Logical, whether to create mutual edges in a directed
 *        ring. It is ignored for undirected graphs.
 * \param circular Logical, if false, the ring will be open (this is
 *        not a real \em ring actually).
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid number of vertices.
 *
 * Time complexity: O(|V|), the
 * number of vertices in the graph.
 *
 * \sa \ref igraph_lattice() for generating more general lattices.
 *
 * \example examples/simple/igraph_ring.c
 */

int igraph_ring(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed,
                igraph_bool_t mutual, igraph_bool_t circular) {

    igraph_vector_t v = IGRAPH_VECTOR_NULL;

    if (n < 0) {
        IGRAPH_ERROR("negative number of vertices", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&v, 1);
    VECTOR(v)[0] = n;

    IGRAPH_CHECK(igraph_lattice(graph, &v, 1, directed, mutual, circular));
    igraph_vector_destroy(&v);

    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * \ingroup generators
 * \function igraph_tree
 * \brief Creates a tree in which almost all vertices have the same number of children.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param n Integer, the number of vertices in the graph.
 * \param children Integer, the number of children of a vertex in the
 *        tree.
 * \param type Constant, gives whether to create a directed tree, and
 *        if this is the case, also its orientation. Possible values:
 *        \clist
 *        \cli IGRAPH_TREE_OUT
 *          directed tree, the edges point
 *          from the parents to their children,
 *        \cli IGRAPH_TREE_IN
 *          directed tree, the edges point from
 *          the children to their parents.
 *        \cli IGRAPH_TREE_UNDIRECTED
 *          undirected tree.
 *        \endclist
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid number of vertices.
 *         \c IGRAPH_INVMODE: invalid mode argument.
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 *
 * \sa \ref igraph_lattice(), \ref igraph_star() for creating other regular
 * structures; \ref igraph_from_prufer() for creating arbitrary trees;
 * \ref igraph_tree_game() for uniform random sampling of trees.
 *
 * \example examples/simple/igraph_tree.c
 */

int igraph_tree(igraph_t *graph, igraph_integer_t n, igraph_integer_t children,
                igraph_tree_mode_t type) {

    igraph_vector_t edges = IGRAPH_VECTOR_NULL;
    long int i, j;
    long int idx = 0;
    long int to = 1;

    if (n < 0 || children <= 0) {
        IGRAPH_ERROR("Invalid number of vertices or children", IGRAPH_EINVAL);
    }
    if (type != IGRAPH_TREE_OUT && type != IGRAPH_TREE_IN &&
        type != IGRAPH_TREE_UNDIRECTED) {
        IGRAPH_ERROR("Invalid mode argument", IGRAPH_EINVMODE);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 2 * (n - 1));

    i = 0;
    if (type == IGRAPH_TREE_OUT) {
        while (idx < 2 * (n - 1)) {
            for (j = 0; j < children && idx < 2 * (n - 1); j++) {
                VECTOR(edges)[idx++] = i;
                VECTOR(edges)[idx++] = to++;
            }
            i++;
        }
    } else {
        while (idx < 2 * (n - 1)) {
            for (j = 0; j < children && idx < 2 * (n - 1); j++) {
                VECTOR(edges)[idx++] = to++;
                VECTOR(edges)[idx++] = i;
            }
            i++;
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, n, type != IGRAPH_TREE_UNDIRECTED));

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * \ingroup generators
 * \function igraph_full
 * \brief Creates a full graph (directed or undirected, with or without loops).
 *
 * </para><para>
 * In a full graph every possible edge is present, every vertex is
 * connected to every other vertex. A full graph in \c igraph should be
 * distinguished from the concept of complete graphs as used in graph theory.
 * If n is a positive integer, then the complete graph K_n on n vertices is
 * the undirected simple graph with the following property. For any distinct
 * pair (u,v) of vertices in K_n, uv (or equivalently vu) is an edge of K_n.
 * In \c igraph, a full graph on n vertices can be K_n, a directed version of
 * K_n, or K_n with at least one loop edge. In any case, if F is a full graph
 * on n vertices as generated by \c igraph, then K_n is a subgraph of the
 * undirected version of F.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param n Integer, the number of vertices in the graph.
 * \param directed Logical, whether to create a directed graph.
 * \param loops Logical, whether to include self-edges (loops).
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid number of vertices.
 *
 * Time complexity: O(|V|+|E|),
 * |V| is the number of vertices,
 * |E| the number of edges in the
 * graph. Of course this is the same as
 * O(|E|)=O(|V||V|)
 * here.
 *
 * \sa \ref igraph_lattice(), \ref igraph_star(), \ref igraph_tree()
 * for creating other regular structures.
 *
 * \example examples/simple/igraph_full.c
 */

int igraph_full(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed,
                igraph_bool_t loops) {

    igraph_vector_t edges = IGRAPH_VECTOR_NULL;
    long int i, j;

    if (n < 0) {
        IGRAPH_ERROR("invalid number of vertices", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

    if (directed && loops) {
        IGRAPH_CHECK(igraph_vector_reserve(&edges, n * n));
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                igraph_vector_push_back(&edges, i); /* reserved */
                igraph_vector_push_back(&edges, j); /* reserved */
            }
        }
    } else if (directed && !loops) {
        IGRAPH_CHECK(igraph_vector_reserve(&edges, n * (n - 1)));
        for (i = 0; i < n; i++) {
            for (j = 0; j < i; j++) {
                igraph_vector_push_back(&edges, i); /* reserved */
                igraph_vector_push_back(&edges, j); /* reserved */
            }
            for (j = i + 1; j < n; j++) {
                igraph_vector_push_back(&edges, i); /* reserved */
                igraph_vector_push_back(&edges, j); /* reserved */
            }
        }
    } else if (!directed && loops) {
        IGRAPH_CHECK(igraph_vector_reserve(&edges, n * (n + 1) / 2));
        for (i = 0; i < n; i++) {
            for (j = i; j < n; j++) {
                igraph_vector_push_back(&edges, i); /* reserved */
                igraph_vector_push_back(&edges, j); /* reserved */
            }
        }
    } else {
        IGRAPH_CHECK(igraph_vector_reserve(&edges, n * (n - 1) / 2));
        for (i = 0; i < n; i++) {
            for (j = i + 1; j < n; j++) {
                igraph_vector_push_back(&edges, i); /* reserved */
                igraph_vector_push_back(&edges, j); /* reserved */
            }
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_full_citation
 * Creates a full citation graph
 *
 * This is a directed graph, where every <code>i->j</code> edge is
 * present if and only if <code>j&lt;i</code>.
 * If the \c directed argument is zero then an undirected graph is
 * created, and it is just a full graph.
 * \param graph Pointer to an uninitialized graph object, the result
 *    is stored here.
 * \param n The number of vertices.
 * \param directed Whether to created a directed graph. If zero an
 *    undirected graph is created.
 * \return Error code.
 *
 * Time complexity: O(|V|^2), as we have many edges.
 */

int igraph_full_citation(igraph_t *graph, igraph_integer_t n,
                         igraph_bool_t directed) {
    igraph_vector_t edges;
    long int i, j, ptr = 0;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, n * (n - 1));
    for (i = 1; i < n; i++) {
        for (j = 0; j < i; j++) {
            VECTOR(edges)[ptr++] = i;
            VECTOR(edges)[ptr++] = j;
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * \function igraph_small
 * \brief Shorthand to create a short graph, giving the edges as arguments.
 *
 * </para><para>
 * This function is handy when a relatively small graph needs to be created.
 * Instead of giving the edges as a vector, they are given simply as
 * arguments and a '-1' needs to be given after the last meaningful
 * edge argument.
 *
 * </para><para>Note that only graphs which have vertices less than
 * the highest value of the 'int' type can be created this way. If you
 * give larger values then the result is undefined.
 *
 * \param graph Pointer to an uninitialized graph object. The result
 *        will be stored here.
 * \param n The number of vertices in the graph; a nonnegative integer.
 * \param directed Logical constant; gives whether the graph should be
 *        directed. Supported values are:
 *        \clist
 *        \cli IGRAPH_DIRECTED
 *          The graph to be created will be \em directed.
 *        \cli IGRAPH_UNDIRECTED
 *          The graph to be created will be \em undirected.
 *        \endclist
 * \param ... The additional arguments giving the edges of the
 *        graph. Don't forget to supply an additional '-1' after the last
 *        (meaningful) argument.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the number
 * of edges in the graph to create.
 *
 * \example examples/simple/igraph_small.c
 */

int igraph_small(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed,
                 ...) {
    igraph_vector_t edges;
    va_list ap;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

    va_start(ap, directed);
    while (1) {
        int num = va_arg(ap, int);
        if (num == -1) {
            break;
        }
        igraph_vector_push_back(&edges, num);
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * \function igraph_extended_chordal_ring
 * Create an extended chordal ring
 *
 * An extended chordal ring is a cycle graph with additional chords
 * connecting its vertices.
 *
 * Each row \c L of the matrix \p W specifies a set of chords to be
 * inserted, in the following way: vertex \c i will connect to a vertex
 * <code>L[(i mod p)]</code> steps ahead of it along the cycle, where
 * \c p is the length of \c L.
 * In other words, vertex \c i will be connected to vertex
 * <code>(i + L[(i mod p)]) mod nodes</code>.
 *
 * </para><para>
 * See also Kotsis, G: Interconnection Topologies for Parallel Processing
 * Systems, PARS Mitteilungen 11, 1-6, 1993.
 *
 * \param graph Pointer to an uninitialized graph object, the result
 *   will be stored here.
 * \param nodes Integer constant, the number of vertices in the
 *   graph. It must be at least 3.
 * \param W The matrix specifying the extra edges. The number of
 *   columns should divide the number of total vertices.
 * \param directed Whether the graph should be directed.
 * \return Error code.
 *
 * \sa \ref igraph_ring(), \ref igraph_lcf(), \ref igraph_lcf_vector()
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the number
 * of edges.
 */

int igraph_extended_chordal_ring(
    igraph_t *graph, igraph_integer_t nodes, const igraph_matrix_t *W,
    igraph_bool_t directed) {
    igraph_vector_t edges;
    long int period = igraph_matrix_ncol(W);
    long int nrow   = igraph_matrix_nrow(W);
    long int i, j, mpos = 0, epos = 0;

    if (nodes < 3) {
        IGRAPH_ERROR("An extended chordal ring has at least 3 nodes", IGRAPH_EINVAL);
    }

    if ((long int)nodes % period != 0) {
        IGRAPH_ERROR("The period (number of columns in W) should divide the "
                     "number of nodes", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 2 * (nodes + nodes * nrow));

    for (i = 0; i < nodes - 1; i++) {
        VECTOR(edges)[epos++] = i;
        VECTOR(edges)[epos++] = i + 1;
    }
    VECTOR(edges)[epos++] = nodes - 1;
    VECTOR(edges)[epos++] = 0;

    if (nrow > 0) {
        for (i = 0; i < nodes; i++) {
            for (j = 0; j < nrow; j++) {
                long int offset = (long int) MATRIX(*W, j, mpos);
                long int v = (i + offset) % nodes;

                if (v < 0) {
                    v += nodes;    /* handle negative offsets */
                }

                VECTOR(edges)[epos++] = i;
                VECTOR(edges)[epos++] = v;

            }
            mpos++; if (mpos == period) {
                mpos = 0;
            }
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_connect_neighborhood
 * \brief Connects every vertex to its neighborhood
 *
 * This function adds new edges to the input graph. Each vertex is connected
 * to all vertices reachable by at most \p order steps from it
 * (unless a connection already existed).  In other words, the \p order power of
 * the graph is computed.
 *
 * </para><para> Note that the input graph is modified in place, no
 * new graph is created. Call \ref igraph_copy() if you want to keep
 * the original graph as well.
 *
 * </para><para> For undirected graphs reachability is always
 * symmetric: if vertex A can be reached from vertex B in at
 * most \p order steps, then the opposite is also true. Only one
 * undirected (A,B) edge will be added in this case.
 * \param graph The input graph, this is the output graph as well.
 * \param order Integer constant, it gives the distance within which
 *    the vertices will be connected to the source vertex.
 * \param mode Constant, it specifies how the neighborhood search is
 *    performed for directed graphs. If \c IGRAPH_OUT then vertices
 *    reachable from the source vertex will be connected, \c IGRAPH_IN
 *    is the opposite. If \c IGRAPH_ALL then the directed graph is
 *    considered as an undirected one.
 * \return Error code.
 *
 * \sa \ref igraph_lattice() uses this function to connect the
 * neighborhood of the vertices.
 *
 * Time complexity: O(|V|*d^k), |V| is the number of vertices in the
 * graph, d is the average degree and k is the \p order argument.
 */

int igraph_connect_neighborhood(igraph_t *graph, igraph_integer_t order,
                                igraph_neimode_t mode) {

    long int no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_t q;
    igraph_vector_t edges;
    long int i, j, in;
    long int *added;
    igraph_vector_t neis;

    if (order < 0) {
        IGRAPH_ERROR("Negative order, cannot connect neighborhood", IGRAPH_EINVAL);
    }

    if (order < 2) {
        IGRAPH_WARNING("Order smaller than two, graph will be unchanged");
    }

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    added = igraph_Calloc(no_of_nodes, long int);
    if (added == 0) {
        IGRAPH_ERROR("Cannot connect neighborhood", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, added);
    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);

    for (i = 0; i < no_of_nodes; i++) {
        added[i] = i + 1;
        igraph_neighbors(graph, &neis, (igraph_integer_t) i, mode);
        in = igraph_vector_size(&neis);
        if (order > 1) {
            for (j = 0; j < in; j++) {
                long int nei = (long int) VECTOR(neis)[j];
                added[nei] = i + 1;
                igraph_dqueue_push(&q, nei);
                igraph_dqueue_push(&q, 1);
            }
        }

        while (!igraph_dqueue_empty(&q)) {
            long int actnode = (long int) igraph_dqueue_pop(&q);
            long int actdist = (long int) igraph_dqueue_pop(&q);
            long int n;
            igraph_neighbors(graph, &neis, (igraph_integer_t) actnode, mode);
            n = igraph_vector_size(&neis);

            if (actdist < order - 1) {
                for (j = 0; j < n; j++) {
                    long int nei = (long int) VECTOR(neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        IGRAPH_CHECK(igraph_dqueue_push(&q, nei));
                        IGRAPH_CHECK(igraph_dqueue_push(&q, actdist + 1));
                        if (mode != IGRAPH_ALL || i < nei) {
                            if (mode == IGRAPH_IN) {
                                IGRAPH_CHECK(igraph_vector_push_back(&edges, nei));
                                IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                            } else {
                                IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                                IGRAPH_CHECK(igraph_vector_push_back(&edges, nei));
                            }
                        }
                    }
                }
            } else {
                for (j = 0; j < n; j++) {
                    long int nei = (long int) VECTOR(neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        if (mode != IGRAPH_ALL || i < nei) {
                            if (mode == IGRAPH_IN) {
                                IGRAPH_CHECK(igraph_vector_push_back(&edges, nei));
                                IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                            } else {
                                IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                                IGRAPH_CHECK(igraph_vector_push_back(&edges, nei));
                            }
                        }
                    }
                }
            }

        } /* while q not empty */
    } /* for i < no_of_nodes */

    igraph_vector_destroy(&neis);
    igraph_dqueue_destroy(&q);
    igraph_free(added);
    IGRAPH_FINALLY_CLEAN(3);

    IGRAPH_CHECK(igraph_add_edges(graph, &edges, 0));

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_de_bruijn
 * \brief Generate a de Bruijn graph.
 *
 * A de Bruijn graph represents relationships between strings. An alphabet
 * of \c m letters are used and strings of length \c n are considered.
 * A vertex corresponds to every possible string and there is a directed edge
 * from vertex \c v to vertex \c w if the string of \c v can be transformed into
 * the string of \c w by removing its first letter and appending a letter to it.
 *
 * </para><para>
 * Please note that the graph will have \c m to the power \c n vertices and
 * even more edges, so probably you don't want to supply too big numbers for
 * \c m and \c n.
 *
 * </para><para>
 * De Bruijn graphs have some interesting properties, please see another source,
 * eg. Wikipedia for details.
 *
 * \param graph Pointer to an uninitialized graph object, the result will be
 *        stored here.
 * \param m Integer, the number of letters in the alphabet.
 * \param n Integer, the length of the strings.
 * \return Error code.
 *
 * \sa \ref igraph_kautz().
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the number of edges.
 */

int igraph_de_bruijn(igraph_t *graph, igraph_integer_t m, igraph_integer_t n) {

    /* m - number of symbols */
    /* n - length of strings */

    long int no_of_nodes, no_of_edges;
    igraph_vector_t edges;
    long int i, j;
    long int mm = m;

    if (m < 0 || n < 0) {
        IGRAPH_ERROR("`m' and `n' should be non-negative in a de Bruijn graph",
                     IGRAPH_EINVAL);
    }

    if (n == 0) {
        return igraph_empty(graph, 1, IGRAPH_DIRECTED);
    }
    if (m == 0) {
        return igraph_empty(graph, 0, IGRAPH_DIRECTED);
    }

    no_of_nodes = (long int) pow(m, n);
    no_of_edges = no_of_nodes * m;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&edges, no_of_edges * 2));

    for (i = 0; i < no_of_nodes; i++) {
        long int basis = (i * mm) % no_of_nodes;
        for (j = 0; j < m; j++) {
            igraph_vector_push_back(&edges, i);
            igraph_vector_push_back(&edges, basis + j);
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, (igraph_integer_t) no_of_nodes,
                               IGRAPH_DIRECTED));

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_kautz
 * \brief Generate a Kautz graph.
 *
 * A Kautz graph is a labeled graph, vertices are labeled by strings
 * of length \c n+1 above an alphabet with \c m+1 letters, with
 * the restriction that every two consecutive letters in the string
 * must be different. There is a directed edge from a vertex \c v to
 * another vertex \c w if it is possible to transform the string of
 * \c v into the string of \c w by removing the first letter and
 * appending a letter to it.
 *
 * </para><para>
 * Kautz graphs have some interesting properties, see eg. Wikipedia
 * for details.
 *
 * </para><para>
 * Vincent Matossian wrote the first version of this function in R,
 * thanks.
 * \param graph Pointer to an uninitialized graph object, the result
 * will be stored here.
 * \param m Integer, \c m+1 is the number of letters in the alphabet.
 * \param n Integer, \c n+1 is the length of the strings.
 * \return Error code.
 *
 * \sa \ref igraph_de_bruijn().
 *
 * Time complexity: O(|V|* [(m+1)/m]^n +|E|), in practice it is more
 * like O(|V|+|E|). |V| is the number of vertices, |E| is the number
 * of edges and \c m and \c n are the corresponding arguments.
 */

int igraph_kautz(igraph_t *graph, igraph_integer_t m, igraph_integer_t n) {

    /* m+1 - number of symbols */
    /* n+1 - length of strings */

    long int mm = m;
    long int no_of_nodes, no_of_edges;
    long int allstrings;
    long int i, j, idx = 0;
    igraph_vector_t edges;
    igraph_vector_long_t digits, table;
    igraph_vector_long_t index1, index2;
    long int actb = 0;
    long int actvalue = 0;

    if (m < 0 || n < 0) {
        IGRAPH_ERROR("`m' and `n' should be non-negative in a Kautz graph",
                     IGRAPH_EINVAL);
    }

    if (n == 0) {
        return igraph_full(graph, m + 1, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    }
    if (m == 0) {
        return igraph_empty(graph, 0, IGRAPH_DIRECTED);
    }

    no_of_nodes = (long int) ((m + 1) * pow(m, n));
    no_of_edges = no_of_nodes * m;
    allstrings = (long int) pow(m + 1, n + 1);

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

    IGRAPH_CHECK(igraph_vector_long_init(&table, n + 1));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &table);
    j = 1;
    for (i = n; i >= 0; i--) {
        VECTOR(table)[i] = j;
        j *= (m + 1);
    }

    IGRAPH_CHECK(igraph_vector_long_init(&digits, n + 1));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &digits);
    IGRAPH_CHECK(igraph_vector_long_init(&index1, (long int) pow(m + 1, n + 1)));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &index1);
    IGRAPH_CHECK(igraph_vector_long_init(&index2, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &index2);

    /* Fill the index tables*/
    while (1) {
        /* at the beginning of the loop, 0:actb contain the valid prefix */
        /* we might need to fill it to get a valid string */
        long int z = 0;
        if (VECTOR(digits)[actb] == 0) {
            z = 1;
        }
        for (actb++; actb <= n; actb++) {
            VECTOR(digits)[actb] = z;
            actvalue += z * VECTOR(table)[actb];
            z = 1 - z;
        }
        actb = n;

        /* ok, we have a valid string now */
        VECTOR(index1)[actvalue] = idx + 1;
        VECTOR(index2)[idx] = actvalue;
        idx++;

        /* finished? */
        if (idx >= no_of_nodes) {
            break;
        }

        /* not yet, we need a valid prefix now */
        while (1) {
            /* try to increase digits at position actb */
            long int next = VECTOR(digits)[actb] + 1;
            if (actb != 0 && VECTOR(digits)[actb - 1] == next) {
                next++;
            }
            if (next <= m) {
                /* ok, no problem */
                actvalue += (next - VECTOR(digits)[actb]) * VECTOR(table)[actb];
                VECTOR(digits)[actb] = next;
                break;
            } else {
                /* bad luck, try the previous digit */
                actvalue -= VECTOR(digits)[actb] * VECTOR(table)[actb];
                actb--;
            }
        }
    }

    IGRAPH_CHECK(igraph_vector_reserve(&edges, no_of_edges * 2));

    /* Now come the edges at last */
    for (i = 0; i < no_of_nodes; i++) {
        long int fromvalue = VECTOR(index2)[i];
        long int lastdigit = fromvalue % (mm + 1);
        long int basis = (fromvalue * (mm + 1)) % allstrings;
        for (j = 0; j <= m; j++) {
            long int tovalue, to;
            if (j == lastdigit) {
                continue;
            }
            tovalue = basis + j;
            to = VECTOR(index1)[tovalue] - 1;
            if (to < 0) {
                continue;
            }
            igraph_vector_push_back(&edges, i);
            igraph_vector_push_back(&edges, to);
        }
    }

    igraph_vector_long_destroy(&index2);
    igraph_vector_long_destroy(&index1);
    igraph_vector_long_destroy(&digits);
    igraph_vector_long_destroy(&table);
    IGRAPH_FINALLY_CLEAN(4);

    IGRAPH_CHECK(igraph_create(graph, &edges, (igraph_integer_t) no_of_nodes,
                               IGRAPH_DIRECTED));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_lcf_vector
 * \brief Create a graph from LCF notation
 *
 * This function is essentially the same as \ref igraph_lcf(), only
 * the way for giving the arguments is different. See \ref
 * igraph_lcf() for details.
 * \param graph Pointer to an uninitialized graph object.
 * \param n Integer constant giving the number of vertices.
 * \param shifts A vector giving the shifts.
 * \param repeats An integer constant giving the number of repeats
 *        for the shifts.
 * \return Error code.
 *
 * \sa \ref igraph_lcf(), \ref igraph_extended_chordal_ring()
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices plus
 * the number of edges.
 */

int igraph_lcf_vector(igraph_t *graph, igraph_integer_t n,
                      const igraph_vector_t *shifts,
                      igraph_integer_t repeats) {

    igraph_vector_t edges;
    long int no_of_shifts = igraph_vector_size(shifts);
    long int ptr = 0, i, sptr = 0;
    long int no_of_nodes = n;
    long int no_of_edges = n + no_of_shifts * repeats;

    if (repeats < 0) {
        IGRAPH_ERROR("number of repeats must be positive", IGRAPH_EINVAL);
    }
    IGRAPH_VECTOR_INIT_FINALLY(&edges, 2 * no_of_edges);

    if (no_of_nodes > 0) {
        /* Create a ring first */
        for (i = 0; i < no_of_nodes; i++) {
            VECTOR(edges)[ptr++] = i;
            VECTOR(edges)[ptr++] = i + 1;
        }
        VECTOR(edges)[ptr - 1] = 0;
    }

    /* Then add the rest */
    while (ptr < 2 * no_of_edges) {
        long int sh = (long int) VECTOR(*shifts)[sptr % no_of_shifts];
        long int from = sptr % no_of_nodes;
        long int to = (no_of_nodes + sptr + sh) % no_of_nodes;
        VECTOR(edges)[ptr++] = from;
        VECTOR(edges)[ptr++] = to;
        sptr++;
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, (igraph_integer_t) no_of_nodes,
                               IGRAPH_UNDIRECTED));
    IGRAPH_CHECK(igraph_simplify(graph, 1 /* true */, 1 /* true */, NULL));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_lcf
 * \brief Create a graph from LCF notation
 *
 * </para><para>
 * LCF is short for Lederberg-Coxeter-Frucht, it is a concise notation for
 * 3-regular Hamiltonian graphs. It consists of three parameters: the
 * number of vertices in the graph, a list of shifts giving additional
 * edges to a cycle backbone, and another integer giving how many times
 * the shifts should be performed. See
 * http://mathworld.wolfram.com/LCFNotation.html for details.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param n Integer, the number of vertices in the graph.
 * \param ... The shifts and the number of repeats for the shifts,
 *        plus an additional 0 to mark the end of the arguments.
 * \return Error code.
 *
 * \sa See \ref igraph_lcf_vector() for a similar function using a
 * vector_t instead of the variable length argument list.
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the number
 * of edges.
 *
 * \example examples/simple/igraph_lcf.c
 */

int igraph_lcf(igraph_t *graph, igraph_integer_t n, ...) {
    igraph_vector_t shifts;
    igraph_integer_t repeats;
    va_list ap;

    IGRAPH_VECTOR_INIT_FINALLY(&shifts, 0);

    va_start(ap, n);
    while (1) {
        int num = va_arg(ap, int);
        if (num == 0) {
            break;
        }
        IGRAPH_CHECK(igraph_vector_push_back(&shifts, num));
    }
    if (igraph_vector_size(&shifts) == 0) {
        repeats = 0;
    } else {
        repeats = (igraph_integer_t) igraph_vector_pop_back(&shifts);
    }

    IGRAPH_CHECK(igraph_lcf_vector(graph, n, &shifts, repeats));
    igraph_vector_destroy(&shifts);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

const igraph_real_t igraph_i_famous_bull[] = {
    5, 5, 0,
    0, 1, 0, 2, 1, 2, 1, 3, 2, 4
};

const igraph_real_t igraph_i_famous_chvatal[] = {
    12, 24, 0,
    5, 6, 6, 7, 7, 8, 8, 9, 5, 9, 4, 5, 4, 8, 2, 8, 2, 6, 0, 6, 0, 9, 3, 9, 3, 7,
    1, 7, 1, 5, 1, 10, 4, 10, 4, 11, 2, 11, 0, 10, 0, 11, 3, 11, 3, 10, 1, 2
};

const igraph_real_t igraph_i_famous_coxeter[] = {
    28, 42, 0,
    0, 1, 0, 2, 0, 7, 1, 4, 1, 13, 2, 3, 2, 8, 3, 6, 3, 9, 4, 5, 4, 12, 5, 6, 5,
    11, 6, 10, 7, 19, 7, 24, 8, 20, 8, 23, 9, 14, 9, 22, 10, 15, 10, 21, 11, 16,
    11, 27, 12, 17, 12, 26, 13, 18, 13, 25, 14, 17, 14, 18, 15, 18, 15, 19, 16, 19,
    16, 20, 17, 20, 21, 23, 21, 26, 22, 24, 22, 27, 23, 25, 24, 26, 25, 27
};

const igraph_real_t igraph_i_famous_cubical[] = {
    8, 12, 0,
    0, 1, 1, 2, 2, 3, 0, 3, 4, 5, 5, 6, 6, 7, 4, 7, 0, 4, 1, 5, 2, 6, 3, 7
};

const igraph_real_t igraph_i_famous_diamond[] = {
    4, 5, 0,
    0, 1, 0, 2, 1, 2, 1, 3, 2, 3
};

const igraph_real_t igraph_i_famous_dodecahedron[] = {
    20, 30, 0,
    0, 1, 0, 4, 0, 5, 1, 2, 1, 6, 2, 3, 2, 7, 3, 4, 3, 8, 4, 9, 5, 10, 5, 11, 6,
    10, 6, 14, 7, 13, 7, 14, 8, 12, 8, 13, 9, 11, 9, 12, 10, 15, 11, 16, 12, 17,
    13, 18, 14, 19, 15, 16, 15, 19, 16, 17, 17, 18, 18, 19
};

const igraph_real_t igraph_i_famous_folkman[] = {
    20, 40, 0,
    0, 5, 0, 8, 0, 10, 0, 13, 1, 7, 1, 9, 1, 12, 1, 14, 2, 6, 2, 8, 2, 11, 2, 13,
    3, 5, 3, 7, 3, 10, 3, 12, 4, 6, 4, 9, 4, 11, 4, 14, 5, 15, 5, 19, 6, 15, 6, 16,
    7, 16, 7, 17, 8, 17, 8, 18, 9, 18, 9, 19, 10, 15, 10, 19, 11, 15, 11, 16, 12,
    16, 12, 17, 13, 17, 13, 18, 14, 18, 14, 19
};

const igraph_real_t igraph_i_famous_franklin[] = {
    12, 18, 0,
    0, 1, 0, 2, 0, 6, 1, 3, 1, 7, 2, 4, 2, 10, 3, 5, 3, 11, 4, 5, 4, 6, 5, 7, 6, 8,
    7, 9, 8, 9, 8, 11, 9, 10, 10, 11
};

const igraph_real_t igraph_i_famous_frucht[] = {
    12, 18, 0,
    0, 1, 0, 2, 0, 11, 1, 3, 1, 6, 2, 5, 2, 10, 3, 4, 3, 6, 4, 8, 4, 11, 5, 9, 5,
    10, 6, 7, 7, 8, 7, 9, 8, 9, 10, 11
};

const igraph_real_t igraph_i_famous_grotzsch[] = {
    11, 20, 0,
    0, 1, 0, 2, 0, 7, 0, 10, 1, 3, 1, 6, 1, 9, 2, 4, 2, 6, 2, 8, 3, 4, 3, 8, 3, 10,
    4, 7, 4, 9, 5, 6, 5, 7, 5, 8, 5, 9, 5, 10
};

const igraph_real_t igraph_i_famous_heawood[] = {
    14, 21, 0,
    0, 1, 0, 5, 0, 13, 1, 2, 1, 10, 2, 3, 2, 7, 3, 4, 3, 12, 4, 5, 4, 9, 5, 6, 6,
    7, 6, 11, 7, 8, 8, 9, 8, 13, 9, 10, 10, 11, 11, 12, 12, 13
};

const igraph_real_t igraph_i_famous_herschel[] = {
    11, 18, 0,
    0, 2, 0, 3, 0, 4, 0, 5, 1, 2, 1, 3, 1, 6, 1, 7, 2, 10, 3, 9, 4, 8, 4, 9, 5, 8,
    5, 10, 6, 8, 6, 9, 7, 8, 7, 10
};

const igraph_real_t igraph_i_famous_house[] = {
    5, 6, 0,
    0, 1, 0, 2, 1, 3, 2, 3, 2, 4, 3, 4
};

const igraph_real_t igraph_i_famous_housex[] = {
    5, 8, 0,
    0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3, 2, 4, 3, 4
};

const igraph_real_t igraph_i_famous_icosahedron[] = {
    12, 30, 0,
    0, 1, 0, 2, 0, 3, 0, 4, 0, 8, 1, 2, 1, 6, 1, 7, 1, 8, 2, 4, 2, 5, 2, 6, 3, 4,
    3, 8, 3, 9, 3, 11, 4, 5, 4, 11, 5, 6, 5, 10, 5, 11, 6, 7, 6, 10, 7, 8, 7, 9, 7,
    10, 8, 9, 9, 10, 9, 11, 10, 11
};

const igraph_real_t igraph_i_famous_krackhardt_kite[] = {
    10, 18, 0,
    0, 1, 0, 2, 0, 3, 0, 5, 1, 3, 1, 4, 1, 6, 2, 3, 2, 5, 3, 4, 3, 5, 3, 6, 4, 6, 5, 6, 5, 7, 6, 7, 7, 8, 8, 9
};

const igraph_real_t igraph_i_famous_levi[] = {
    30, 45, 0,
    0, 1, 0, 7, 0, 29, 1, 2, 1, 24, 2, 3, 2, 11, 3, 4, 3, 16, 4, 5, 4, 21, 5, 6, 5,
    26, 6, 7, 6, 13, 7, 8, 8, 9, 8, 17, 9, 10, 9, 22, 10, 11, 10, 27, 11, 12, 12,
    13, 12, 19, 13, 14, 14, 15, 14, 23, 15, 16, 15, 28, 16, 17, 17, 18, 18, 19, 18,
    25, 19, 20, 20, 21, 20, 29, 21, 22, 22, 23, 23, 24, 24, 25, 25, 26, 26, 27, 27,
    28, 28, 29
};

const igraph_real_t igraph_i_famous_mcgee[] = {
    24, 36, 0,
    0, 1, 0, 7, 0, 23, 1, 2, 1, 18, 2, 3, 2, 14, 3, 4, 3, 10, 4, 5, 4, 21, 5, 6, 5,
    17, 6, 7, 6, 13, 7, 8, 8, 9, 8, 20, 9, 10, 9, 16, 10, 11, 11, 12, 11, 23, 12,
    13, 12, 19, 13, 14, 14, 15, 15, 16, 15, 22, 16, 17, 17, 18, 18, 19, 19, 20, 20,
    21, 21, 22, 22, 23
};

const igraph_real_t igraph_i_famous_meredith[] = {
    70, 140, 0,
    0, 4, 0, 5, 0, 6, 1, 4, 1, 5, 1, 6, 2, 4, 2, 5, 2, 6, 3, 4, 3, 5, 3, 6, 7, 11,
    7, 12, 7, 13, 8, 11, 8, 12, 8, 13, 9, 11, 9, 12, 9, 13, 10, 11, 10, 12, 10, 13,
    14, 18, 14, 19, 14, 20, 15, 18, 15, 19, 15, 20, 16, 18, 16, 19, 16, 20, 17, 18,
    17, 19, 17, 20, 21, 25, 21, 26, 21, 27, 22, 25, 22, 26, 22, 27, 23, 25, 23, 26,
    23, 27, 24, 25, 24, 26, 24, 27, 28, 32, 28, 33, 28, 34, 29, 32, 29, 33, 29, 34,
    30, 32, 30, 33, 30, 34, 31, 32, 31, 33, 31, 34, 35, 39, 35, 40, 35, 41, 36, 39,
    36, 40, 36, 41, 37, 39, 37, 40, 37, 41, 38, 39, 38, 40, 38, 41, 42, 46, 42, 47,
    42, 48, 43, 46, 43, 47, 43, 48, 44, 46, 44, 47, 44, 48, 45, 46, 45, 47, 45, 48,
    49, 53, 49, 54, 49, 55, 50, 53, 50, 54, 50, 55, 51, 53, 51, 54, 51, 55, 52, 53,
    52, 54, 52, 55, 56, 60, 56, 61, 56, 62, 57, 60, 57, 61, 57, 62, 58, 60, 58, 61,
    58, 62, 59, 60, 59, 61, 59, 62, 63, 67, 63, 68, 63, 69, 64, 67, 64, 68, 64, 69,
    65, 67, 65, 68, 65, 69, 66, 67, 66, 68, 66, 69, 2, 50, 1, 51, 9, 57, 8, 58, 16,
    64, 15, 65, 23, 36, 22, 37, 30, 43, 29, 44, 3, 21, 7, 24, 14, 31, 0, 17, 10,
    28, 38, 42, 35, 66, 59, 63, 52, 56, 45, 49
};

const igraph_real_t igraph_i_famous_noperfectmatching[] = {
    16, 27, 0,
    0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3, 2, 4, 3, 4, 4, 5, 5, 6, 5, 7, 6, 12, 6, 13,
    7, 8, 7, 9, 8, 9, 8, 10, 8, 11, 9, 10, 9, 11, 10, 11, 12, 13, 12, 14, 12, 15,
    13, 14, 13, 15, 14, 15
};

const igraph_real_t igraph_i_famous_nonline[] = {
    50, 72, 0,
    0, 1, 0, 2, 0, 3, 4, 6, 4, 7, 5, 6, 5, 7, 6, 7, 7, 8, 9, 11, 9, 12, 9, 13, 10,
    11, 10, 12, 10, 13, 11, 12, 11, 13, 12, 13, 14, 15, 15, 16, 15, 17, 16, 17, 16,
    18, 17, 18, 18, 19, 20, 21, 20, 22, 20, 23, 21, 22, 21, 23, 21, 24, 22, 23, 22,
    24, 24, 25, 26, 27, 26, 28, 26, 29, 27, 28, 27, 29, 27, 30, 27, 31, 28, 29, 28,
    30, 28, 31, 30, 31, 32, 34, 32, 35, 32, 36, 33, 34, 33, 35, 33, 37, 34, 35, 36,
    37, 38, 39, 38, 40, 38, 43, 39, 40, 39, 41, 39, 42, 39, 43, 40, 41, 41, 42, 42,
    43, 44, 45, 44, 46, 45, 46, 45, 47, 46, 47, 46, 48, 47, 48, 47, 49, 48, 49
};

const igraph_real_t igraph_i_famous_octahedron[] = {
    6, 12, 0,
    0, 1, 0, 2, 1, 2, 3, 4, 3, 5, 4, 5, 0, 3, 0, 5, 1, 3, 1, 4, 2, 4, 2, 5
};

const igraph_real_t igraph_i_famous_petersen[] = {
    10, 15, 0,
    0, 1, 0, 4, 0, 5, 1, 2, 1, 6, 2, 3, 2, 7, 3, 4, 3, 8, 4, 9, 5, 7, 5, 8, 6, 8, 6, 9, 7, 9
};

const igraph_real_t igraph_i_famous_robertson[] = {
    19, 38, 0,
    0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12,
    12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 0, 18, 0, 4, 4, 9, 9, 13, 13,
    17, 2, 17, 2, 6, 6, 10, 10, 15, 0, 15, 1, 8, 8, 16, 5, 16, 5, 12, 1, 12, 7, 18,
    7, 14, 3, 14, 3, 11, 11, 18
};

const igraph_real_t igraph_i_famous_smallestcyclicgroup[] = {
    9, 15, 0,
    0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 1, 2, 1, 3, 1, 7, 1, 8, 2, 5, 2, 6, 2, 7, 3, 8,
    4, 5, 6, 7
};

const igraph_real_t igraph_i_famous_tetrahedron[] = {
    4, 6, 0,
    0, 3, 1, 3, 2, 3, 0, 1, 1, 2, 0, 2
};

const igraph_real_t igraph_i_famous_thomassen[] = {
    34, 52, 0,
    0, 2, 0, 3, 1, 3, 1, 4, 2, 4, 5, 7, 5, 8, 6, 8, 6, 9, 7, 9, 10, 12, 10, 13, 11,
    13, 11, 14, 12, 14, 15, 17, 15, 18, 16, 18, 16, 19, 17, 19, 9, 19, 4, 14, 24,
    25, 25, 26, 20, 26, 20, 21, 21, 22, 22, 23, 23, 27, 27, 28, 28, 29, 29, 30, 30,
    31, 31, 32, 32, 33, 24, 33, 5, 24, 6, 25, 7, 26, 8, 20, 0, 20, 1, 21, 2, 22, 3,
    23, 10, 27, 11, 28, 12, 29, 13, 30, 15, 30, 16, 31, 17, 32, 18, 33
};

const igraph_real_t igraph_i_famous_tutte[] = {
    46, 69, 0,
    0, 10, 0, 11, 0, 12, 1, 2, 1, 7, 1, 19, 2, 3, 2, 41, 3, 4, 3, 27, 4, 5, 4, 33,
    5, 6, 5, 45, 6, 9, 6, 29, 7, 8, 7, 21, 8, 9, 8, 22, 9, 24, 10, 13, 10, 14, 11,
    26, 11, 28, 12, 30, 12, 31, 13, 15, 13, 21, 14, 15, 14, 18, 15, 16, 16, 17, 16,
    20, 17, 18, 17, 23, 18, 24, 19, 25, 19, 40, 20, 21, 20, 22, 22, 23, 23, 24, 25,
    26, 25, 38, 26, 34, 27, 28, 27, 39, 28, 34, 29, 30, 29, 44, 30, 35, 31, 32, 31,
    35, 32, 33, 32, 42, 33, 43, 34, 36, 35, 37, 36, 38, 36, 39, 37, 42, 37, 44, 38,
    40, 39, 41, 40, 41, 42, 43, 43, 45, 44, 45
};

const igraph_real_t igraph_i_famous_uniquely3colorable[] = {
    12, 22, 0,
    0, 1, 0, 3, 0, 6, 0, 8, 1, 4, 1, 7, 1, 9, 2, 3, 2, 6, 2, 7, 2, 9, 2, 11, 3, 4,
    3, 10, 4, 5, 4, 11, 5, 6, 5, 7, 5, 8, 5, 10, 8, 11, 9, 10
};

const igraph_real_t igraph_i_famous_walther[] = {
    25, 31, 0,
    0, 1, 1, 2, 1, 8, 2, 3, 2, 13, 3, 4, 3, 16, 4, 5, 5, 6, 5, 19, 6, 7, 6, 20, 7,
    21, 8, 9, 8, 13, 9, 10, 9, 22, 10, 11, 10, 20, 11, 12, 13, 14, 14, 15, 14, 23,
    15, 16, 15, 17, 17, 18, 18, 19, 18, 24, 20, 24, 22, 23, 23, 24
};

const igraph_real_t igraph_i_famous_zachary[] = {
    34, 78, 0,
    0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8,
    0, 10, 0, 11, 0, 12, 0, 13, 0, 17, 0, 19, 0, 21, 0, 31,
    1, 2, 1, 3, 1, 7, 1, 13, 1, 17, 1, 19, 1, 21, 1, 30,
    2, 3, 2, 7, 2, 27, 2, 28, 2, 32, 2, 9, 2, 8, 2, 13,
    3, 7, 3, 12, 3, 13, 4, 6, 4, 10, 5, 6, 5, 10, 5, 16,
    6, 16, 8, 30, 8, 32, 8, 33, 9, 33, 13, 33, 14, 32, 14, 33,
    15, 32, 15, 33, 18, 32, 18, 33, 19, 33, 20, 32, 20, 33,
    22, 32, 22, 33, 23, 25, 23, 27, 23, 32, 23, 33, 23, 29,
    24, 25, 24, 27, 24, 31, 25, 31, 26, 29, 26, 33, 27, 33,
    28, 31, 28, 33, 29, 32, 29, 33, 30, 32, 30, 33, 31, 32, 31, 33,
    32, 33
};

static int igraph_i_famous(igraph_t *graph, const igraph_real_t *data) {
    long int no_of_nodes = (long int) data[0];
    long int no_of_edges = (long int) data[1];
    igraph_bool_t directed = (igraph_bool_t) data[2];
    igraph_vector_t edges;

    igraph_vector_view(&edges, data + 3, 2 * no_of_edges);
    IGRAPH_CHECK(igraph_create(graph, &edges, (igraph_integer_t) no_of_nodes,
                               directed));
    return 0;
}

/**
 * \function igraph_famous
 * \brief Create a famous graph by simply providing its name
 *
 * </para><para>
 * The name of the graph can be simply supplied as a string.
 * Note that this function creates graphs which don't take any parameters,
 * there are separate functions for graphs with parameters, eg. \ref
 * igraph_full() for creating a full graph.
 *
 * </para><para>
 * The following graphs are supported:
 * \clist
 *   \cli Bull
 *           The bull graph, 5 vertices, 5 edges, resembles the
 *           head of a bull if drawn properly.
 *   \cli Chvatal
 *           This is the smallest triangle-free graph that is
 *           both 4-chromatic and 4-regular. According to the Grunbaum
 *           conjecture there exists an m-regular, m-chromatic graph
 *           with n vertices for every m>1 and n>2. The Chvatal graph
 *           is an example for m=4 and n=12. It has 24 edges.
 *   \cli Coxeter
 *           A non-Hamiltonian cubic symmetric graph with 28
 *           vertices and 42 edges.
 *   \cli Cubical
 *           The Platonic graph of the cube. A convex regular
 *           polyhedron with 8 vertices and 12 edges.
 *   \cli Diamond
 *           A graph with 4 vertices and 5 edges, resembles a
 *           schematic diamond if drawn properly.
 *   \cli Dodecahedral, Dodecahedron
 *           Another Platonic solid
 *           with 20 vertices and 30 edges.
 *   \cli Folkman
 *           The semisymmetric graph with minimum number of
 *           vertices, 20 and 40 edges. A semisymmetric graph is
 *           regular, edge transitive and not vertex transitive.
 *   \cli Franklin
 *           This is a graph whose embedding to the Klein
 *           bottle can be colored with six colors, it is a
 *           counterexample to the necessity of the Heawood
 *           conjecture on a Klein bottle. It has 12 vertices and 18
 *           edges.
 *   \cli Frucht
 *           The Frucht Graph is the smallest cubical graph
 *           whose automorphism group consists only of the identity
 *           element. It has 12 vertices and 18 edges.
 *   \cli Grotzsch
 *           The Grötzsch graph is a triangle-free graph with
 *           11 vertices, 20 edges, and chromatic number 4. It is named after
 *           German mathematician Herbert Grötzsch, and its existence
 *           demonstrates that the assumption of planarity is necessary in
 *           Grötzsch's theorem that every triangle-free planar
 *           graph is 3-colorable.
 *   \cli Heawood
 *           The Heawood graph is an undirected graph with 14
 *           vertices and 21 edges. The graph is cubic, and all cycles in the
 *           graph have six or more edges. Every smaller cubic graph has shorter
 *           cycles, so this graph is the 6-cage, the smallest cubic graph of
 *           girth 6.
 *   \cli Herschel
 *           The Herschel graph is the smallest
 *           nonhamiltonian polyhedral graph. It is the
 *           unique such graph on 11 nodes, and has 18 edges.
 *   \cli House
 *           The house graph is a 5-vertex, 6-edge graph, the
 *           schematic draw of a house if drawn properly, basically a
 *           triangle on top of a square.
 *   \cli HouseX
 *           The same as the house graph with an X in the square. 5
 *           vertices and 8 edges.
 *   \cli Icosahedral, Icosahedron
 *           A Platonic solid with 12
 *           vertices and 30 edges.
 *   \cli Krackhardt_Kite
 *           A social network with 10 vertices and 18 edges.
 *           Krackhardt, D. Assessing the Political Landscape:
 *           Structure, Cognition, and Power in Organizations.
 *           Admin. Sci. Quart. 35, 342-369, 1990.
 *   \cli Levi
 *           The graph is a 4-arc transitive cubic graph, it has
 *           30 vertices and 45 edges.
 *   \cli McGee
 *           The McGee graph is the unique 3-regular 7-cage
 *           graph, it has 24 vertices and 36 edges.
 *   \cli Meredith
 *           The Meredith graph is a quartic graph on 70
 *           nodes and 140 edges that is a counterexample to the conjecture that
 *           every 4-regular 4-connected graph is Hamiltonian.
 *   \cli Noperfectmatching
 *           A connected graph with 16 vertices and
 *           27 edges containing no perfect matching. A matching in a graph
 *           is a set of pairwise non-incident edges; that is, no two edges
 *           share a common vertex. A perfect matching is a matching
 *           which covers all vertices of the graph.
 *   \cli Nonline
 *           A graph whose connected components are the 9
 *           graphs whose presence as a vertex-induced subgraph in a
 *           graph makes a nonline graph. It has 50 vertices and 72 edges.
 *   \cli Octahedral, Octahedron
 *           Platonic solid with 6
 *           vertices and 12 edges.
 *   \cli Petersen
 *           A 3-regular graph with 10 vertices and 15 edges. It is
 *           the smallest hypohamiltonian graph, ie. it is
 *           non-hamiltonian but removing any single vertex from it makes it
 *           Hamiltonian.
 *   \cli Robertson
 *           The unique (4,5)-cage graph, ie. a 4-regular
 *           graph of girth 5. It has 19 vertices and 38 edges.
 *   \cli Smallestcyclicgroup
 *           A smallest nontrivial graph
 *           whose automorphism group is cyclic. It has 9 vertices and
 *           15 edges.
 *   \cli Tetrahedral, Tetrahedron
 *           Platonic solid with 4
 *           vertices and 6 edges.
 *   \cli Thomassen
 *           The smallest hypotraceable graph,
 *           on 34 vertices and 52 edges. A hypotracable graph does
 *           not contain a Hamiltonian path but after removing any
 *           single vertex from it the remainder always contains a
 *           Hamiltonian path. A graph containing a Hamiltonian path
 *           is called traceable.
 *   \cli Tutte
 *           Tait's Hamiltonian graph conjecture states that
 *           every 3-connected 3-regular planar graph is Hamiltonian.
 *           This graph is a counterexample. It has 46 vertices and 69
 *           edges.
 *   \cli Uniquely3colorable
 *           Returns a 12-vertex, triangle-free
 *           graph with chromatic number 3 that is uniquely
 *           3-colorable.
 *   \cli Walther
 *           An identity graph with 25 vertices and 31
 *           edges. An identity graph has a single graph automorphism,
 *           the trivial one.
 *   \cli Zachary
 *           Social network of friendships between 34 members of a
 *           karate club at a US university in the 1970s. See
 *           W. W. Zachary, An information flow model for conflict and
 *           fission in small groups, Journal of Anthropological
 *           Research 33, 452-473 (1977).
 * \endclist
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param name Character constant, the name of the graph to be
 *     created, it is case insensitive.
 * \return Error code, IGRAPH_EINVAL if there is no graph with the
 *     given name.
 *
 * \sa Other functions for creating graph structures:
 * \ref igraph_ring(), \ref igraph_tree(), \ref igraph_lattice(), \ref
 * igraph_full().
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the number
 * of edges in the graph.
 */

int igraph_famous(igraph_t *graph, const char *name) {

    if (!strcasecmp(name, "bull")) {
        return igraph_i_famous(graph, igraph_i_famous_bull);
    } else if (!strcasecmp(name, "chvatal")) {
        return igraph_i_famous(graph, igraph_i_famous_chvatal);
    } else if (!strcasecmp(name, "coxeter")) {
        return igraph_i_famous(graph, igraph_i_famous_coxeter);
    } else if (!strcasecmp(name, "cubical")) {
        return igraph_i_famous(graph, igraph_i_famous_cubical);
    } else if (!strcasecmp(name, "diamond")) {
        return igraph_i_famous(graph, igraph_i_famous_diamond);
    } else if (!strcasecmp(name, "dodecahedral") ||
               !strcasecmp(name, "dodecahedron")) {
        return igraph_i_famous(graph, igraph_i_famous_dodecahedron);
    } else if (!strcasecmp(name, "folkman")) {
        return igraph_i_famous(graph, igraph_i_famous_folkman);
    } else if (!strcasecmp(name, "franklin")) {
        return igraph_i_famous(graph, igraph_i_famous_franklin);
    } else if (!strcasecmp(name, "frucht")) {
        return igraph_i_famous(graph, igraph_i_famous_frucht);
    } else if (!strcasecmp(name, "grotzsch")) {
        return igraph_i_famous(graph, igraph_i_famous_grotzsch);
    } else if (!strcasecmp(name, "heawood")) {
        return igraph_i_famous(graph, igraph_i_famous_heawood);
    } else if (!strcasecmp(name, "herschel")) {
        return igraph_i_famous(graph, igraph_i_famous_herschel);
    } else if (!strcasecmp(name, "house")) {
        return igraph_i_famous(graph, igraph_i_famous_house);
    } else if (!strcasecmp(name, "housex")) {
        return igraph_i_famous(graph, igraph_i_famous_housex);
    } else if (!strcasecmp(name, "icosahedral") ||
               !strcasecmp(name, "icosahedron")) {
        return igraph_i_famous(graph, igraph_i_famous_icosahedron);
    } else if (!strcasecmp(name, "krackhardt_kite")) {
        return igraph_i_famous(graph, igraph_i_famous_krackhardt_kite);
    } else if (!strcasecmp(name, "levi")) {
        return igraph_i_famous(graph, igraph_i_famous_levi);
    } else if (!strcasecmp(name, "mcgee")) {
        return igraph_i_famous(graph, igraph_i_famous_mcgee);
    } else if (!strcasecmp(name, "meredith")) {
        return igraph_i_famous(graph, igraph_i_famous_meredith);
    } else if (!strcasecmp(name, "noperfectmatching")) {
        return igraph_i_famous(graph, igraph_i_famous_noperfectmatching);
    } else if (!strcasecmp(name, "nonline")) {
        return igraph_i_famous(graph, igraph_i_famous_nonline);
    } else if (!strcasecmp(name, "octahedral") ||
               !strcasecmp(name, "octahedron")) {
        return igraph_i_famous(graph, igraph_i_famous_octahedron);
    } else if (!strcasecmp(name, "petersen")) {
        return igraph_i_famous(graph, igraph_i_famous_petersen);
    } else if (!strcasecmp(name, "robertson")) {
        return igraph_i_famous(graph, igraph_i_famous_robertson);
    } else if (!strcasecmp(name, "smallestcyclicgroup")) {
        return igraph_i_famous(graph, igraph_i_famous_smallestcyclicgroup);
    } else if (!strcasecmp(name, "tetrahedral") ||
               !strcasecmp(name, "tetrahedron")) {
        return igraph_i_famous(graph, igraph_i_famous_tetrahedron);
    } else if (!strcasecmp(name, "thomassen")) {
        return igraph_i_famous(graph, igraph_i_famous_thomassen);
    } else if (!strcasecmp(name, "tutte")) {
        return igraph_i_famous(graph, igraph_i_famous_tutte);
    } else if (!strcasecmp(name, "uniquely3colorable")) {
        return igraph_i_famous(graph, igraph_i_famous_uniquely3colorable);
    } else if (!strcasecmp(name, "walther")) {
        return igraph_i_famous(graph, igraph_i_famous_walther);
    } else if (!strcasecmp(name, "zachary")) {
        return igraph_i_famous(graph, igraph_i_famous_zachary);
    } else {
        IGRAPH_ERROR("Unknown graph, see documentation", IGRAPH_EINVAL);
    }

    return 0;
}

/**
 * \function igraph_adjlist
 * Create a graph from an adjacency list
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


/**
 * \ingroup generators
 * \function igraph_from_prufer
 * \brief Generates a tree from a Pr&uuml;fer sequence
 *
 * A Pr&uuml;fer sequence is a unique sequence of integers associated
 * with a labelled tree. A tree on n vertices can be represented by a
 * sequence of n-2 integers, each between 0 and n-1 (inclusive).
 *
 * The algorithm used by this function is based on
 * Paulius Micikevi&ccaron;ius, Saverio Caminiti, Narsingh Deo:
 * Linear-time Algorithms for Encoding Trees as Sequences of Node Labels
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param prufer The Pr&uuml;fer sequence
 * \return Error code:
 *          \clist
 *          \cli IGRAPH_ENOMEM
 *             there is not enough memory to perform the operation.
 *          \cli IGRAPH_EINVAL
 *             invalid Pr&uuml;fer sequence given
 *          \endclist
 *
 * \sa \ref igraph_tree(), \ref igraph_tree_game()
 *
 */

int igraph_from_prufer(igraph_t *graph, const igraph_vector_int_t *prufer) {
    igraph_vector_int_t degree;
    igraph_vector_t edges;
    long n;
    long i, k;
    long u, v; /* vertices */
    long ec;

    n = igraph_vector_int_size(prufer) + 2;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&degree, n); /* initializes vector to zeros */
    IGRAPH_VECTOR_INIT_FINALLY(&edges, 2 * (n - 1));

    /* build out-degree vector (i.e. number of child vertices) and verify Prufer sequence */
    for (i = 0; i < n - 2; ++i) {
        long u = VECTOR(*prufer)[i];
        if (u >= n || u < 0) {
            IGRAPH_ERROR("Invalid Prufer sequence", IGRAPH_EINVAL);
        }
        VECTOR(degree)[u] += 1;
    }

    v = 0;  /* initialize v now, in case Prufer sequence is empty */
    k = 0;  /* index into the Prufer vector */
    ec = 0; /* index into the edges vector */
    for (i = 0; i < n; ++i) {
        u = i;

        while (k < n - 2 && u <= i && (VECTOR(degree)[u] == 0)) {
            /* u is a leaf here */

            v = VECTOR(*prufer)[k]; /* parent of u */

            /* add edge */
            VECTOR(edges)[ec++] = v;
            VECTOR(edges)[ec++] = u;

            k += 1;

            VECTOR(degree)[v] -= 1;

            u = v;
        }

        if (k == n - 2) {
            break;
        }
    }

    /* find u for last edge, v is already set */
    for (u = i + 1; u < n; ++u)
        if ((VECTOR(degree)[u] == 0) && u != v) {
            break;
        }

    /* add last edge */
    VECTOR(edges)[ec++] = v;
    VECTOR(edges)[ec++] = u;

    IGRAPH_CHECK(igraph_create(graph, &edges, (igraph_integer_t) n, /* directed = */ 0));

    igraph_vector_destroy(&edges);
    igraph_vector_int_destroy(&degree);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}
