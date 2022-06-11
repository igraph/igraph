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

#include "math/safe_intop.h"

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
 * \sa \ref igraph_square_lattice(), \ref igraph_star(), \ref igraph_kary_tree()
 * for creating other regular structures.
 *
 * \example examples/simple/igraph_full.c
 */
igraph_error_t igraph_full(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed,
                igraph_bool_t loops) {

    igraph_vector_int_t edges = IGRAPH_VECTOR_NULL;
    igraph_integer_t no_of_edges2;
    igraph_integer_t i, j;

    if (n < 0) {
        IGRAPH_ERROR("Invalid number of vertices.", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

    if (directed && loops) {
        /* ecount = n * n */
        IGRAPH_SAFE_MULT(n, n, &no_of_edges2);
        IGRAPH_SAFE_MULT(no_of_edges2, 2, &no_of_edges2);
        IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_of_edges2));
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                igraph_vector_int_push_back(&edges, i); /* reserved */
                igraph_vector_int_push_back(&edges, j); /* reserved */
            }
        }
    } else if (directed && !loops) {
        /* ecount = n * (n - 1) */
        IGRAPH_SAFE_MULT(n, n - 1, &no_of_edges2);
        IGRAPH_SAFE_MULT(no_of_edges2, 2, &no_of_edges2);
        IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_of_edges2));
        for (i = 0; i < n; i++) {
            for (j = 0; j < i; j++) {
                igraph_vector_int_push_back(&edges, i); /* reserved */
                igraph_vector_int_push_back(&edges, j); /* reserved */
            }
            for (j = i + 1; j < n; j++) {
                igraph_vector_int_push_back(&edges, i); /* reserved */
                igraph_vector_int_push_back(&edges, j); /* reserved */
            }
        }
    } else if (!directed && loops) {
        /* ecount = n * (n + 1) / 2 */
        IGRAPH_SAFE_ADD(n, 1, &no_of_edges2);
        IGRAPH_SAFE_MULT(n, no_of_edges2, &no_of_edges2);
        IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_of_edges2));
        for (i = 0; i < n; i++) {
            for (j = i; j < n; j++) {
                igraph_vector_int_push_back(&edges, i); /* reserved */
                igraph_vector_int_push_back(&edges, j); /* reserved */
            }
        }
    } else {
        /* ecount = n * (n - 1) / 2 */
        IGRAPH_SAFE_MULT(n, n - 1, &no_of_edges2);
        IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_of_edges2));
        for (i = 0; i < n; i++) {
            for (j = i + 1; j < n; j++) {
                igraph_vector_int_push_back(&edges, i); /* reserved */
                igraph_vector_int_push_back(&edges, j); /* reserved */
            }
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_full_multipartite
 * \brief Create a full multipartite network.
 *
 * A multipartite network contains two or more types and connections
 * are only possible between two vertices in different types.
 *
 * \param graph Pointer to an igraph_t object, the graph will be
 *   created here.
 * \param types Pointer to an int vector. If not a null pointer,
 *   it contains information about the vertex types
 * \param n Pointer to an int vector, the number of types
 * \param directed Boolean, whether to create a directed graph.
 * \param mode A constant that gives the type of connections for
 *   directed graphs. If \c IGRAPH_OUT, then edges point from vertices
 *   of low-index vertices to high-index vertices; if \c
 *   IGRAPH_IN, then the opposite direction is realized; if \c
 *   IGRAPH_ALL, then mutual edges will be created.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 *
 * \sa \ref igraph_full_bipartite() for bipartite full graphs.
 */
igraph_error_t igraph_full_multipartite(igraph_t *graph,
                          igraph_vector_int_t *types,
                          const igraph_vector_int_t *n,
                          igraph_bool_t directed,
                          igraph_neimode_t mode) {
    
    igraph_integer_t no_of_edges = 0;
    igraph_integer_t no_of_nodes;
    igraph_vector_int_t edges;
    igraph_vector_int_t n_acc;
    igraph_integer_t ptr = 0;

    igraph_integer_t nn = igraph_vector_int_size(n);

    if (nn == 0) {
        igraph_empty(graph, 0, directed);
        if (types) {
            igraph_vector_int_clear(types);
        }
        return IGRAPH_SUCCESS;
    }

    if (nn == 1) {
        igraph_integer_t num = VECTOR(*n)[0];
        igraph_empty(graph, num, directed);
        if (types) {
            IGRAPH_CHECK(igraph_vector_int_resize(types, num));
            igraph_vector_int_null(types);
            for (igraph_integer_t i = 0; i < num; i++) {
                VECTOR(*types)[0] = 0;
            }
        }
        return IGRAPH_SUCCESS;
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&n_acc, nn);
    VECTOR(n_acc)[0] = 0;
    for (igraph_integer_t i = 1; i < nn; i++) {
        IGRAPH_SAFE_ADD(VECTOR(n_acc)[i-1], VECTOR(*n)[i-1],
                    &VECTOR(n_acc)[i]);
    }
    IGRAPH_SAFE_ADD(VECTOR(n_acc)[nn-1], VECTOR(*n)[nn-1], &no_of_nodes);

    for (igraph_integer_t i = 0; i < nn; i++) {
        igraph_integer_t v = VECTOR(*n)[i];
        igraph_integer_t partial_sum;
        IGRAPH_SAFE_ADD(no_of_nodes, -v, &partial_sum);
        IGRAPH_SAFE_MULT(partial_sum, v, &partial_sum);
        IGRAPH_SAFE_ADD(no_of_edges, partial_sum, &no_of_edges);
    }

    if (directed && mode == IGRAPH_ALL) {
        IGRAPH_SAFE_MULT(no_of_edges, 2, &no_of_edges);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges);

    for (igraph_integer_t from = 0; from < nn-1; from++) {
        igraph_integer_t edge_from = VECTOR(n_acc)[from];
        for (igraph_integer_t i = 0; i < VECTOR(*n)[from]; i++) {
            for (igraph_integer_t to = from+1; to < nn; to++) {
                igraph_integer_t edge_to = VECTOR(n_acc)[to];
                for (igraph_integer_t j = 0; j < VECTOR(*n)[to]; j++) {
                    if (!directed || mode == IGRAPH_OUT) {
                        VECTOR(edges)[ptr++] = edge_from;
                        VECTOR(edges)[ptr++] = edge_to++;
                    } else if (mode == IGRAPH_IN) {
                        VECTOR(edges)[ptr++] = edge_to++;
                        VECTOR(edges)[ptr++] = edge_from;
                    } else {
                        VECTOR(edges)[ptr++] = edge_from;
                        VECTOR(edges)[ptr++] = edge_to;
                        VECTOR(edges)[ptr++] = edge_to++;
                        VECTOR(edges)[ptr++] = edge_from;
                    }
                }
            }
            edge_from++;
        }
    }
    
    IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, directed));

    if (types) {
        IGRAPH_CHECK(igraph_vector_int_resize(types, no_of_nodes));
        igraph_vector_int_null(types);
        if (no_of_nodes > 0) {
            igraph_integer_t v = 1;
            VECTOR(*types)[0] = 0;
            for (igraph_integer_t i = 1; i < no_of_nodes; i++) {
                if (v < nn && i == VECTOR(n_acc)[v]) {
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
 * \function igraph_full_citation
 * \brief Creates a full citation graph.
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
igraph_error_t igraph_full_citation(igraph_t *graph, igraph_integer_t n,
                         igraph_bool_t directed) {
    igraph_vector_int_t edges;
    igraph_integer_t i, j, ptr = 0;

    if (n < 0) {
        IGRAPH_ERROR("Invalid number of vertices.", IGRAPH_EINVAL);
    }

    {
        igraph_integer_t no_of_edges2;
        IGRAPH_SAFE_MULT(n, n-1, &no_of_edges2);
        IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges2);
    }

    for (i = 1; i < n; i++) {
        for (j = 0; j < i; j++) {
            VECTOR(edges)[ptr++] = i;
            VECTOR(edges)[ptr++] = j;
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}
