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
#include "igraph_memory.h"
#include "igraph_operators.h"

#include "core/interruption.h"

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
 * \brief Arbitrary dimensional square lattices.
 *
 * Creates d-dimensional square lattices of the given size. Optionally,
 * the lattice can be made periodic, and the neighbors within a given
 * graph distance can be connected.
 *
 * </para><para>
 * In the zero-dimensional case, the singleton graph is returned.
 *
 * </para><para>
 * The vertices of the resulting graph are ordered such that the
 * index of the vertex at position <code>(i_0, i_1, i_2, ..., i_d)</code>
 * in a lattice of size <code>(n_0, n_1, ..., n_d)</code> will be
 * <code>i_0 + n_0 * i_1 + n_0 * n_1 * i_2 + ...</code>.
 *
 * \param graph An uninitialized graph object.
 * \param dimvector Vector giving the sizes of the lattice in each of
 *        its dimensions. The dimension of the lattice will be the
 *        same as the length of this vector.
 * \param nei Integer value giving the distance (number of steps)
 *        within which two vertices will be connected.
 * \param directed Boolean, whether to create a directed graph. 
 *        If the \c mutual and \c circular arguments are not set to true,
 *        edges will be directed from lower-index vertices towards
 *        higher-index ones.
 * \param mutual Boolean, if the graph is directed this gives whether
 *        to create all connections as mutual.
 * \param circular Boolean, defines whether the generated lattice is
 *        periodic.
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid (negative)
 *         dimension vector.
 *
 * Time complexity: If \p nei is less than two then it is O(|V|+|E|) (as
 * far as I remember), |V| and |E| are the number of vertices
 * and edges in the generated graph. Otherwise it is O(|V|*d^k+|E|), d
 * is the average degree of the graph, k is the \p nei argument.
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
        IGRAPH_ERROR("Lattice creation failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, coords);
    weights = igraph_Calloc(dims, long int);
    if (weights == 0) {
        IGRAPH_ERROR("Lattice creation failed", IGRAPH_ENOMEM);
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
