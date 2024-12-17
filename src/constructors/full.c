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

#include "igraph_interface.h"

#include "core/interruption.h"
#include "math/safe_intop.h"

/**
 * \ingroup generators
 * \function igraph_full
 * \brief Creates a full graph (complete graph).
 *
 * In a full graph every possible edge is present: every vertex is
 * connected to every other vertex. \a igraph generalizes the usual
 * concept of complete graphs in graph theory to graphs with self-loops
 * as well as to directed graphs.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param n Integer, the number of vertices in the graph.
 * \param directed Logical, whether to create a directed graph.
 * \param loops Logical, whether to include self-edges (loops).
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid number of vertices.
 *
 * Time complexity: O(|V|^2) = O(|E|),
 * where |V| is the number of vertices and |E| is the number of edges.
 *
 * \sa \ref igraph_square_lattice(), \ref igraph_star(), \ref igraph_kary_tree()
 * for creating other regular structures.
 *
 * \example examples/simple/igraph_full.c
 */
igraph_error_t igraph_full(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed,
                igraph_bool_t loops) {

    igraph_vector_int_t edges = IGRAPH_VECTOR_NULL;
    igraph_integer_t no_of_edges2;

    if (n < 0) {
        IGRAPH_ERROR("Invalid number of vertices.", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

    if (directed && loops) {
        /* ecount = n * n */
        IGRAPH_SAFE_MULT(n, n, &no_of_edges2);
        IGRAPH_SAFE_MULT(no_of_edges2, 2, &no_of_edges2);
        IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_of_edges2));
        for (igraph_integer_t i = 0; i < n; i++) {
            for (igraph_integer_t j = 0; j < n; j++) {
                igraph_vector_int_push_back(&edges, i); /* reserved */
                igraph_vector_int_push_back(&edges, j); /* reserved */
            }
            IGRAPH_ALLOW_INTERRUPTION();
        }
    } else if (directed && !loops) {
        /* ecount = n * (n - 1) */
        IGRAPH_SAFE_MULT(n, n - 1, &no_of_edges2);
        IGRAPH_SAFE_MULT(no_of_edges2, 2, &no_of_edges2);
        IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_of_edges2));
        for (igraph_integer_t i = 0; i < n; i++) {
            for (igraph_integer_t j = 0; j < i; j++) {
                igraph_vector_int_push_back(&edges, i); /* reserved */
                igraph_vector_int_push_back(&edges, j); /* reserved */
            }
            for (igraph_integer_t j = i + 1; j < n; j++) {
                igraph_vector_int_push_back(&edges, i); /* reserved */
                igraph_vector_int_push_back(&edges, j); /* reserved */
            }
            IGRAPH_ALLOW_INTERRUPTION();
        }
    } else if (!directed && loops) {
        /* ecount = n * (n + 1) / 2 */
        IGRAPH_SAFE_ADD(n, 1, &no_of_edges2);
        IGRAPH_SAFE_MULT(n, no_of_edges2, &no_of_edges2);
        IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_of_edges2));
        for (igraph_integer_t i = 0; i < n; i++) {
            for (igraph_integer_t j = i; j < n; j++) {
                igraph_vector_int_push_back(&edges, i); /* reserved */
                igraph_vector_int_push_back(&edges, j); /* reserved */
            }
            IGRAPH_ALLOW_INTERRUPTION();
        }
    } else {
        /* ecount = n * (n - 1) / 2 */
        IGRAPH_SAFE_MULT(n, n - 1, &no_of_edges2);
        IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_of_edges2));
        for (igraph_integer_t i = 0; i < n; i++) {
            for (igraph_integer_t j = i + 1; j < n; j++) {
                igraph_vector_int_push_back(&edges, i); /* reserved */
                igraph_vector_int_push_back(&edges, j); /* reserved */
            }
            IGRAPH_ALLOW_INTERRUPTION();
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_full_multipartite
 * \brief Creates a full multipartite graph.
 *
 * A multipartite graph contains two or more types of vertices and connections
 * are only possible between two vertices of different types. This function
 * creates a complete multipartite graph.
 *
 * \param graph Pointer to an uninitialized graph object, the graph will be
 *   created here.
 * \param types Pointer to an integer vector. If not a null pointer,
 *   the type of each vertex will be stored here.
 * \param n Pointer to an integer vector, the number of vertices
 *   of each type.
 * \param directed Boolean, whether to create a directed graph.
 * \param mode A constant that gives the type of connections for
 *   directed graphs. If \c IGRAPH_OUT, then edges point from vertices
 *   of low-index vertices to high-index vertices; if
 *   \c IGRAPH_IN, then the opposite direction is realized;
 *   \c IGRAPH_ALL, then mutual edges will be created.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 *
 * \sa \ref igraph_full_bipartite() for complete bipartite graphs,
 * \ref igraph_turan() for Turán graphs.
 */
igraph_error_t igraph_full_multipartite(igraph_t *graph,
                          igraph_vector_int_t *types,
                          const igraph_vector_int_t *n,
                          igraph_bool_t directed,
                          igraph_neimode_t mode) {

    igraph_vector_int_t edges;
    igraph_vector_int_t n_acc;

    igraph_integer_t no_of_types = igraph_vector_int_size(n);

    if (no_of_types == 0) {
        IGRAPH_CHECK(igraph_empty(graph, 0, directed));
        if (types) {
            igraph_vector_int_clear(types);
        }
        return IGRAPH_SUCCESS;
    }

    if (igraph_vector_int_min(n) < 0) {
        IGRAPH_ERROR("Number of vertices must not be negative in any partition.", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&n_acc, no_of_types+1);
    VECTOR(n_acc)[0] = 0;
    for (igraph_integer_t i = 1; i < no_of_types+1; i++) {
        IGRAPH_SAFE_ADD(VECTOR(n_acc)[i-1], VECTOR(*n)[i-1],
                    &VECTOR(n_acc)[i]);
    }

    igraph_integer_t no_of_edges2 = 0;

    for (igraph_integer_t i = 0; i < no_of_types; i++) {
        igraph_integer_t v = VECTOR(*n)[i];
        igraph_integer_t partial_sum = VECTOR(n_acc)[no_of_types] - v;
        IGRAPH_SAFE_MULT(partial_sum, v, &partial_sum);
        IGRAPH_SAFE_ADD(no_of_edges2, partial_sum, &no_of_edges2);
    }

    if (directed && mode == IGRAPH_ALL) {
        IGRAPH_SAFE_MULT(no_of_edges2, 2, &no_of_edges2);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges2);

    igraph_integer_t ptr = 0;

    for (igraph_integer_t from_type = 0; from_type < no_of_types-1; from_type++) {
        igraph_integer_t edge_from = VECTOR(n_acc)[from_type];
        for (igraph_integer_t i = 0; i < VECTOR(*n)[from_type]; i++) {
            for (igraph_integer_t to_type = from_type+1; to_type < no_of_types; to_type++) {
                igraph_integer_t edge_to = VECTOR(n_acc)[to_type];
                for (igraph_integer_t j = 0; j < VECTOR(*n)[to_type]; j++) {
                    if (!directed || mode == IGRAPH_OUT) {
                        VECTOR(edges)[ptr++] = edge_from;
                        VECTOR(edges)[ptr++] = edge_to;
                    } else if (mode == IGRAPH_IN) {
                        VECTOR(edges)[ptr++] = edge_to;
                        VECTOR(edges)[ptr++] = edge_from;
                    } else {
                        VECTOR(edges)[ptr++] = edge_from;
                        VECTOR(edges)[ptr++] = edge_to;
                        VECTOR(edges)[ptr++] = edge_to;
                        VECTOR(edges)[ptr++] = edge_from;
                    }
                    edge_to++;
                }
            }
            edge_from++;
            IGRAPH_ALLOW_INTERRUPTION();
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, VECTOR(n_acc)[no_of_types], directed));

    if (types) {
        IGRAPH_CHECK(igraph_vector_int_resize(types, VECTOR(n_acc)[no_of_types]));
        if (VECTOR(n_acc)[no_of_types] > 0) {
            igraph_integer_t v = 1;
            for (igraph_integer_t i = 0; i < VECTOR(n_acc)[no_of_types]; i++) {
                if (i == VECTOR(n_acc)[v]) {
                    v++;
                }
                VECTOR(*types)[i] = v-1;
            }
        }
    }

    igraph_vector_int_destroy(&edges);
    igraph_vector_int_destroy(&n_acc);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_turan
 * \brief Creates a Turán graph.
 *
 * Turán graphs are complete multipartite graphs with the property
 * that the sizes of the partitions are as close to equal as possible.
 *
 * </para><para>
 * The Turán graph with \p n vertices and \p r partitions is the densest
 * graph on \p n  vertices that does not contain a clique of size
 * <code>r+1</code>.
 *
 * </para><para>
 * This function generates undirected graphs. The null graph is
 * returned when the number of vertices is zero. A complete graph is
 * returned if the number of partitions is greater than the number of
 * vertices.
 *
 * \param graph Pointer to an igraph_t object, the graph will be
 *   created here.
 * \param types Pointer to an integer vector. If not a null pointer,
 *   the type (partition index) of each vertex will be stored here.
 * \param n Integer, the number of vertices in the graph.
 * \param r Integer, the number of partitions of the graph, must be
 *   positive.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 *
 * \sa \ref igraph_full_multipartite() for full multipartite graphs.
 */
igraph_error_t igraph_turan(igraph_t *graph,
                            igraph_vector_int_t *types,
                            igraph_integer_t n,
                            igraph_integer_t r) {
    igraph_integer_t quotient;
    igraph_integer_t remainder;
    igraph_vector_int_t subsets;

    if (n < 0) {
        IGRAPH_ERRORF("Number of vertices must not be negative, got %" IGRAPH_PRId ".", IGRAPH_EINVAL, n);
    }

    if (r <= 0) {
        IGRAPH_ERRORF("Number of partitions must be positive, got %" IGRAPH_PRId ".", IGRAPH_EINVAL, r);
    }

    if (n == 0) {
        IGRAPH_CHECK(igraph_empty(graph, 0, IGRAPH_UNDIRECTED));
        if (types) {
            igraph_vector_int_clear(types);
        }
        return IGRAPH_SUCCESS;
    }

    if (r > n) {
        r = n;
    }

    quotient = n / r;
    remainder = n % r;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&subsets, r);

    igraph_vector_int_fill(&subsets, quotient);
    for (igraph_integer_t i = 0; i < remainder; i++) {
        VECTOR(subsets)[i]++;
    }

    IGRAPH_CHECK(igraph_full_multipartite(graph, types, &subsets,
            IGRAPH_UNDIRECTED, IGRAPH_ALL));
    igraph_vector_int_destroy(&subsets);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_full_citation
 * \brief Creates a full citation graph (a complete directed acyclic graph).
 *
 * This is a directed graph, where every <code>i->j</code> edge is
 * present if and only if <code>j&lt;i</code>.
 * If the \p directed argument is false then an undirected graph is
 * created, and it is just a complete graph.
 *
 * \param graph Pointer to an uninitialized graph object, the result
 *    is stored here.
 * \param n The number of vertices.
 * \param directed Whether to created a directed graph. If false an
 *    undirected graph is created.
 * \return Error code.
 *
 * \sa \ref igraph_full()
 *
 * Time complexity: O(|V|^2) = O(|E|),
 * where |V| is the number of vertices and |E| is the number of edges.
 */
igraph_error_t igraph_full_citation(igraph_t *graph, igraph_integer_t n,
                         igraph_bool_t directed) {
    igraph_vector_int_t edges;
    igraph_integer_t ptr = 0;

    if (n < 0) {
        IGRAPH_ERROR("Invalid number of vertices.", IGRAPH_EINVAL);
    }

    {
        igraph_integer_t no_of_edges2;
        IGRAPH_SAFE_MULT(n, n-1, &no_of_edges2);
        IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges2);
    }

    for (igraph_integer_t i = 1; i < n; i++) {
        for (igraph_integer_t j = 0; j < i; j++) {
            VECTOR(edges)[ptr++] = i;
            VECTOR(edges)[ptr++] = j;
        }
        IGRAPH_ALLOW_INTERRUPTION();
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
