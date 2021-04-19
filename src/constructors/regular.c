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

    coords = IGRAPH_CALLOC(dims, long int);
    if (coords == 0) {
        IGRAPH_ERROR("Lattice creation failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, coords);
    weights = IGRAPH_CALLOC(dims, long int);
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
    IGRAPH_FREE(coords);
    IGRAPH_FREE(weights);
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
 * \function igraph_extended_chordal_ring
 * \brief Create an extended chordal ring.
 *
 * An extended chordal ring is a cycle graph with additional chords
 * connecting its vertices.
 *
 * Each row \c L of the matrix \p W specifies a set of chords to be
 * inserted, in the following way: vertex \c i will connect to a vertex
 * <code>L[(i mod p)]</code> steps ahead of it along the cycle, where
 * \c p is the length of \c L.
 * In other words, vertex \c i will be connected to vertex
 * <code>(i + L[(i mod p)]) mod nodes</code>. If multiple edges are
 * defined in this way, this will output a non-simple graph. The result
 * can be simplified using \ref igraph_simplify().
 *
 * </para><para>
 * See also Kotsis, G: Interconnection Topologies for Parallel Processing
 * Systems, PARS Mitteilungen 11, 1-6, 1993. The igraph extended chordal
 * rings are not identical to the ones in the paper. In igraph
 * the matrix specifies which edges to add. In the paper, a condition is
 * specified which should simultaneously hold between two endpoints and
 * the reverse endpoints.
 *
 * \param graph Pointer to an uninitialized graph object, the result
 *   will be stored here.
 * \param nodes Integer constant, the number of vertices in the
 *   graph. It must be at least 3.
 * \param W The matrix specifying the extra edges. The number of
 *   columns should divide the number of total vertices. The elements
 *   are allowed to be negative.
 * \param directed Whether the graph should be directed.
 * \return Error code.
 *
 * \sa \ref igraph_ring(), \ref igraph_lcf(), \ref igraph_lcf_vector().
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
