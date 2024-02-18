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
#include "math/safe_intop.h"

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
 * \sa \ref igraph_square_lattice(), \ref igraph_ring(), \ref igraph_kary_tree()
 * for creating other regular structures.
 *
 * \example examples/simple/igraph_star.c
 */
igraph_error_t igraph_star(igraph_t *graph, igraph_integer_t n, igraph_star_mode_t mode,
                igraph_integer_t center) {

    igraph_vector_int_t edges = IGRAPH_VECTOR_NULL;
    igraph_integer_t i;

    if (n < 0) {
        IGRAPH_ERROR("Invalid number of vertices.", IGRAPH_EINVVID);
    }
    if (center < 0 || center > n - 1) {
        IGRAPH_ERROR("Invalid center vertex.", IGRAPH_EINVAL);
    }
    if (mode != IGRAPH_STAR_OUT && mode != IGRAPH_STAR_IN &&
        mode != IGRAPH_STAR_MUTUAL && mode != IGRAPH_STAR_UNDIRECTED) {
        IGRAPH_ERROR("Invalid star mode.", IGRAPH_EINVMODE);
    }

    if (mode != IGRAPH_STAR_MUTUAL) {
        igraph_integer_t no_of_edges2;
        IGRAPH_SAFE_MULT(n-1, 2, &no_of_edges2);
        IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges2);
    } else {
        igraph_integer_t no_of_edges2;
        IGRAPH_SAFE_MULT(n-1, 4, &no_of_edges2);
        IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges2);
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
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup generators
 * \function igraph_wheel
 * \brief Creates a \em wheel graph, a union of a star and a cycle graph.
 *
 * A wheel graph on \p n vertices can be thought of as a wheel with
 * <code>n - 1</code> spokes. The cycle graph part makes up the rim,
 * while the star graph part adds the spokes.
 *
 * </para><para>
 * Note that the two and three-vertex wheel graphs are non-simple:
 * The two-vertex wheel graph contains a self-loop, while the three-vertex
 * wheel graph contains parallel edges (a 1-cycle and a 2-cycle, respectively).
 *
 * \param graph Pointer to an uninitialized graph object, this will
 *        be the result.
 * \param n Integer constant, the number of vertices in the graph.
 * \param mode Constant, gives the type of the star graph to
 *        create. Possible values:
 *        \clist
 *        \cli IGRAPH_WHEEL_OUT
 *          directed wheel graph, edges point
 *          \em from the center to the other vertices.
 *        \cli IGRAPH_WHEEL_IN
 *          directed wheel graph, edges point
 *          \em to the center from the other vertices.
 *        \cli IGRAPH_WHEEL_MUTUAL
 *          directed wheel graph with mutual edges.
 *        \cli IGRAPH_WHEEL_UNDIRECTED
 *          an undirected wheel graph is
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
 * \sa \ref igraph_square_lattice(), \ref igraph_ring(), \ref igraph_star(),
 * \ref igraph_kary_tree() for creating other regular structures.
 *
 */

igraph_error_t igraph_wheel(igraph_t *graph, igraph_integer_t n, igraph_wheel_mode_t mode,
                igraph_integer_t center) {

    igraph_star_mode_t star_mode;
    igraph_vector_int_t rim_edges = IGRAPH_VECTOR_NULL;
    igraph_integer_t i;

    /* Firstly creates a star by the function \ref igraph_star() and makes
     * use of its existing input parameter checking ability, it can check
     * "Invalid number of vertices" and "Invalid center vertex". */
    switch (mode)
    {
        case IGRAPH_WHEEL_OUT:
            star_mode = IGRAPH_STAR_OUT;
            break;
        case IGRAPH_WHEEL_IN:
            star_mode = IGRAPH_STAR_IN;
            break;
        case IGRAPH_WHEEL_MUTUAL:
            star_mode = IGRAPH_STAR_MUTUAL;
            break;
        case IGRAPH_WHEEL_UNDIRECTED:
            star_mode = IGRAPH_STAR_UNDIRECTED;
            break;
        default:
            IGRAPH_ERROR("Invalid wheel graph mode.", IGRAPH_EINVMODE);
    }

    IGRAPH_CHECK(igraph_star(graph, n, star_mode, center));

    /* If n <= 1, wheel graph is identical with star graph,
     * no further processing is needed. */
    if (n <= 1) {
        return IGRAPH_SUCCESS;
    }

    /* Register the star for deallocation in case of error flow before
     * the entire wheel is successfully created. */
    IGRAPH_FINALLY(igraph_destroy, graph);

    /* Add edges to the rim. As the rim (or cycle) has n - 1 vertices,
     * it will have n - 1 edges. For MUTUAL mode, number of edges
     * will be double. */
    if (mode == IGRAPH_WHEEL_MUTUAL) {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&rim_edges, 4 * (n-1));
    } else {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&rim_edges, 2 * (n-1));
    }

    /* Assign first n-1 edges (MUTUAL will be handled later). */
    for (i = 0; i < n-2; i++) {
        if ( i < center ) {
            VECTOR(rim_edges)[2 * i] = i;
            if ( i + 1 < center ) {
                VECTOR(rim_edges)[2 * i + 1] = i + 1;
            } else {
                VECTOR(rim_edges)[2 * i + 1] = i + 2;
            }
        } else {
            VECTOR(rim_edges)[2 * i] = i + 1;
            VECTOR(rim_edges)[2 * i + 1] = i + 2;
        }
    }

    /* Assign the last edge (MUTUAL will be handled later). */
    if ( n - 2 < center ) {
        VECTOR(rim_edges)[2 * n - 4] = n - 2;
    } else {
        VECTOR(rim_edges)[2 * n - 4] = n - 1;
    }
    if ( center > 0 ) {
        VECTOR(rim_edges)[2 * n - 3] = 0;
    } else {
        VECTOR(rim_edges)[2 * n - 3] = 1;
    }

    /* For MUTUAL mode, add reverse-direction edges. */
    if (mode == IGRAPH_WHEEL_MUTUAL) {
        for (i=0; i < 2 * (n-1); i++) {
            VECTOR(rim_edges)[4 * (n-1) - 1 - i] = VECTOR(rim_edges)[i];
        }
    }

    /* Combine the rim into the star to make it a wheel graph. */
    IGRAPH_CHECK(igraph_add_edges(graph, &rim_edges, NULL));

    igraph_vector_int_destroy(&rim_edges);

    /* 2 instead of 1 because the star graph is registered before. */
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup generators
 * \function igraph_lattice
 * \brief Arbitrary dimensional square lattices (deprecated).
 *
 * \deprecated-by igraph_square_lattice 0.10.0
 */
igraph_error_t igraph_lattice(igraph_t *graph, const igraph_vector_int_t *dimvector,
                   igraph_integer_t nei, igraph_bool_t directed, igraph_bool_t mutual,
                   igraph_bool_t circular) {
    igraph_vector_bool_t periodic;

    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&periodic, igraph_vector_int_size(dimvector));
    igraph_vector_bool_fill(&periodic, circular);

    IGRAPH_CHECK(igraph_square_lattice(graph, dimvector, nei, directed, mutual, &periodic));

    igraph_vector_bool_destroy(&periodic);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup generators
 * \function igraph_square_lattice
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
 * index of the vertex at position <code>(i_1, i_2, i_3, ..., i_d)</code>
 * in a lattice of size <code>(n_1, n_2, ..., n_d)</code> will be
 * <code>i_1 + n_1 * i_2 + n_1 * n_2 * i_3 + ...</code>.
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
 * \param periodic Boolean vector, defines whether the generated lattice is
 *        periodic along each dimension. The length of this vector must match
 *        the length of \p dimvector. This parameter may also be \c NULL, which
 *        implies that the lattice will not be periodic.
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid (negative) dimension vector or mismatch
 *         between the length of the dimension vector and the periodicity vector.
 *
 * Time complexity: If \p nei is less than two then it is O(|V|+|E|) (as
 * far as I remember), |V| and |E| are the number of vertices
 * and edges in the generated graph. Otherwise it is O(|V|*d^k+|E|), d
 * is the average degree of the graph, k is the \p nei argument.
 */
igraph_error_t igraph_square_lattice(
    igraph_t *graph, const igraph_vector_int_t *dimvector, igraph_integer_t nei,
    igraph_bool_t directed, igraph_bool_t mutual, const igraph_vector_bool_t *periodic
) {

    igraph_integer_t dims = igraph_vector_int_size(dimvector);
    igraph_integer_t no_of_nodes;
    igraph_vector_int_t edges = IGRAPH_VECTOR_NULL;
    igraph_integer_t *coords, *weights;
    igraph_integer_t i, j;
    int carry, pos;
    int iter = 0;

    if (igraph_vector_int_any_smaller(dimvector, 0)) {
        IGRAPH_ERROR("Invalid dimension vector.", IGRAPH_EINVAL);
    }

    if (periodic && igraph_vector_bool_size(periodic) != dims) {
        IGRAPH_ERRORF(
            "Length of periodicity vector must match the length of the "
            "dimension vector (%" IGRAPH_PRId ").",
            IGRAPH_EINVAL, dims
        );
    }

    /* compute no. of nodes in overflow-safe manner */
    IGRAPH_CHECK(igraph_i_safe_vector_int_prod(dimvector, &no_of_nodes));

    /* init coords & weights */

    coords = IGRAPH_CALLOC(dims, igraph_integer_t);
    IGRAPH_CHECK_OOM(coords, "Lattice creation failed.");
    IGRAPH_FINALLY(igraph_free, coords);

    weights = IGRAPH_CALLOC(dims, igraph_integer_t);
    IGRAPH_CHECK_OOM(weights, "Lattice creation failed.");
    IGRAPH_FINALLY(igraph_free, weights);

    if (dims > 0) {
        weights[0] = 1;
        for (i = 1; i < dims; i++) {
            weights[i] = weights[i - 1] * VECTOR(*dimvector)[i - 1];
        }
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    if (mutual && directed) {
        igraph_integer_t no_of_edges2;
        IGRAPH_SAFE_MULT(no_of_nodes, dims, &no_of_edges2);
        IGRAPH_SAFE_MULT(no_of_edges2, 2, &no_of_edges2);
        IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_of_edges2));
    } else {
        igraph_integer_t no_of_edges2;
        IGRAPH_SAFE_MULT(no_of_nodes, dims, &no_of_edges2);
        IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_of_edges2));
    }

#define IS_PERIODIC(dim) ((periodic && VECTOR(*periodic)[dim]))

    for (i = 0; i < no_of_nodes; i++) {
        IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 10);

        /* Connect the current node to the "next" node along each dimension */
        for (j = 0; j < dims; j++) {
            igraph_bool_t is_periodic = IS_PERIODIC(j);

            if (is_periodic|| coords[j] != VECTOR(*dimvector)[j] - 1) {
                igraph_integer_t new_nei;
                if (coords[j] != VECTOR(*dimvector)[j] - 1) {
                    new_nei = i + weights[j] + 1;
                } else {
                    new_nei = i - (VECTOR(*dimvector)[j] - 1) * weights[j] + 1;
                }
                if (new_nei != i + 1 &&
                    (VECTOR(*dimvector)[j] != 2 || coords[j] != 1 || directed)) {
                    igraph_vector_int_push_back(&edges, i); /* reserved */
                    igraph_vector_int_push_back(&edges, new_nei - 1); /* reserved */
                }
            } /* if is_periodic || coords[j] */
            if (mutual && directed && (is_periodic || coords[j] != 0)) {
                igraph_integer_t new_nei;
                if (coords[j] != 0) {
                    new_nei = i - weights[j] + 1;
                } else {
                    new_nei = i + (VECTOR(*dimvector)[j] - 1) * weights[j] + 1;
                }
                if (new_nei != i + 1 &&
                    (VECTOR(*dimvector)[j] != 2 || !is_periodic)) {
                    igraph_vector_int_push_back(&edges, i); /* reserved */
                    igraph_vector_int_push_back(&edges, new_nei - 1); /* reserved */
                }
            } /* if is_periodic || coords[0] */
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

    IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, directed));
    if (nei >= 2) {
        IGRAPH_CHECK(igraph_connect_neighborhood(graph, nei, IGRAPH_ALL));
    }

    /* clean up */
    IGRAPH_FREE(coords);
    IGRAPH_FREE(weights);
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup generators
 * \function igraph_ring
 * \brief Creates a \em cycle graph or a \em path graph.
 *
 * A circular ring on \c n vertices is commonly known in graph
 * theory as the cycle graph, and often denoted by <code>C_n</code>.
 * Removing a single edge from the cycle graph <code>C_n</code> results
 * in the path graph <code>P_n</code>. This function can generate both.
 *
 * </para><para>
 * When \p n is 1 or 2, the result may not be a simple graph:
 * the one-cycle contains a self-loop and the undirected or reciprocally
 * connected directed two-cycle contains parallel edges.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param n The number of vertices in the graph.
 * \param directed Logical, whether to create a directed graph.
 *        All edges will be oriented in the same direction along
 *        the cycle or path.
 * \param mutual Logical, whether to create mutual edges in directed
 *        graphs. It is ignored for undirected graphs.
 * \param circular Logical, whether to create a closed ring (a cycle)
 *        or an open path.
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid number of vertices.
 *
 * Time complexity: O(|V|), the number of vertices in the graph.
 *
 * \sa \ref igraph_lattice() for generating more general lattices.
 *
 * \example examples/simple/igraph_ring.c
 */
igraph_error_t igraph_ring(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed,
                igraph_bool_t mutual, igraph_bool_t circular) {

    igraph_vector_int_t edges;
    igraph_integer_t no_of_edges, no_of_edges2;
    igraph_integer_t i;

    if (n < 0) {
        IGRAPH_ERRORF("The number of vertices must be non-negative, got %" IGRAPH_PRId ".", IGRAPH_EINVAL, n);
    }

    if (n == 0) {
        return igraph_empty(graph, 0, directed);
    }

    no_of_edges = circular ? n : n-1;
    if (directed && mutual) {
        IGRAPH_SAFE_MULT(no_of_edges, 2, &no_of_edges);
    }
    IGRAPH_SAFE_MULT(no_of_edges, 2, &no_of_edges2);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges2);

    if (directed && mutual) {
        for (i=0; i < n-1; ++i) {
            VECTOR(edges)[4*i]   = i;
            VECTOR(edges)[4*i+1] = i+1;
            VECTOR(edges)[4*i+2] = i+1;
            VECTOR(edges)[4*i+3] = i;
        }
        if (circular) {
            /* Now i == n-1 */
            VECTOR(edges)[4*i]   = i;
            VECTOR(edges)[4*i+1] = 0;
            VECTOR(edges)[4*i+2] = 0;
            VECTOR(edges)[4*i+3] = i;
        }
    } else {
        for (i=0; i < n-1; ++i) {
            VECTOR(edges)[2*i]   = i;
            VECTOR(edges)[2*i+1] = i+1;
        }
        if (circular) {
            /* Now i == n-1 */
            VECTOR(edges)[2*i]   = i;
            VECTOR(edges)[2*i+1] = 0;
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup generators
 * \function igraph_kary_tree
 * \brief Creates a k-ary tree in which almost all vertices have k children.
 *
 * To obtain a completely symmetric tree with \c l layers, where each
 * vertex has precisely \p children descendants, use
 * <code>n = (children^(l+1) - 1) / (children - 1)</code>.
 * Such trees are often called <code>k</code>-ary trees, where \c k refers
 * to the number of children.
 *
 * </para><para>
 * Note that for <code>n=0</code>, the null graph is returned,
 * which is not considered to be a tree by \ref igraph_is_tree().
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
 *          from the parents to their children.
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
 * \example examples/simple/igraph_kary_tree.c
 */
igraph_error_t igraph_kary_tree(igraph_t *graph, igraph_integer_t n, igraph_integer_t children,
                igraph_tree_mode_t type) {

    igraph_vector_int_t edges = IGRAPH_VECTOR_NULL;
    igraph_integer_t i, j;
    igraph_integer_t idx = 0;
    igraph_integer_t to = 1;

    if (n < 0) {
        IGRAPH_ERROR("Number of vertices cannot be negative.", IGRAPH_EINVAL);
    }
    if (children <= 0) {
        IGRAPH_ERROR("Number of children must be positive.", IGRAPH_EINVAL);
    }
    if (type != IGRAPH_TREE_OUT && type != IGRAPH_TREE_IN &&
        type != IGRAPH_TREE_UNDIRECTED) {
        IGRAPH_ERROR("Invalid tree orientation type.", IGRAPH_EINVMODE);
    }

    {
        igraph_integer_t no_of_edges2;
        if (n > 0) {
            IGRAPH_SAFE_MULT(n-1, 2, &no_of_edges2);
        } else {
            no_of_edges2 = 0;
        }
        IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges2);
    }

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

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup generators
 * \function igraph_tree
 * \brief Creates a k-ary tree in which almost all vertices have k children (deprecated alias).
 *
 * \deprecated-by igraph_kary_tree 0.10.0
 */
igraph_error_t igraph_tree(igraph_t *graph, igraph_integer_t n, igraph_integer_t children,
                igraph_tree_mode_t type) {
    return igraph_kary_tree(graph, n, children, type);
}

/**
 * \ingroup generators
 * \function igraph_symmetric_tree
 * \brief Creates a symmetric tree with the specified number of branches at each level.
 *
 * This function creates a tree in which all vertices at distance \c d from the
 * root have \p branching_counts[d] children.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param branches Vector detailing the number of branches at each level.
 * \param type Constant, gives whether to create a directed tree, and
 *        if this is the case, also its orientation. Possible values:
 *        \clist
 *        \cli IGRAPH_TREE_OUT
 *          directed tree, the edges point
 *          from the parents to their children.
 *        \cli IGRAPH_TREE_IN
 *          directed tree, the edges point from
 *          the children to their parents.
 *        \cli IGRAPH_TREE_UNDIRECTED
 *          undirected tree.
 *        \endclist
 * \return Error code:
 *         \c IGRAPH_INVMODE: invalid mode argument.
 *         \c IGRAPH_EINVAL: invalid number of children.
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 *
 * \sa \ref igraph_kary_tree(), \ref igraph_regular_tree() and \ref igraph_star()
 * for creating other regular tree structures;
 * \ref igraph_from_prufer() for creating arbitrary trees;
 * \ref igraph_tree_game() for uniform random sampling of trees.
 *
 * \example examples/simple/igraph_symmetric_tree.c
 */

igraph_error_t igraph_symmetric_tree(igraph_t *graph, const igraph_vector_int_t *branches,
                igraph_tree_mode_t type) {

    igraph_vector_int_t edges;
    igraph_integer_t j, k, temp, no_of_nodes, idx, parent, child, level_end;
    igraph_integer_t branching_counts_size = igraph_vector_int_size(branches);

    if (type != IGRAPH_TREE_OUT && type != IGRAPH_TREE_IN && type != IGRAPH_TREE_UNDIRECTED) {
        IGRAPH_ERROR("Invalid tree orientation type.", IGRAPH_EINVMODE);
    }
    if (!igraph_vector_int_empty(branches) && igraph_vector_int_min(branches) <= 0) {
        IGRAPH_ERROR("The number of branches must be positive at each level.", IGRAPH_EINVAL);
    }

    /* Compute the number of vertices in the tree. */
    no_of_nodes = 1;
    temp = 1;
    for (j = 0; j < branching_counts_size; ++j) {
        IGRAPH_SAFE_MULT(temp, VECTOR(*branches)[j], &temp);
        IGRAPH_SAFE_ADD(no_of_nodes, temp, &no_of_nodes);
    }

    /* Trees have precisely |E| = |V| - 1 edges. */
    {
        igraph_integer_t no_of_edges2;
        IGRAPH_SAFE_MULT(no_of_nodes - 1, 2, &no_of_edges2);
        IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges2);
    }

    idx = 0;

    /* Current parent and child vertex ids.
     * parent -> child edges will be added. */
    child = 1;
    parent = 0;
    for (k = 0; k < branching_counts_size; ++k) {
        level_end = child; /* points to one past the last vertex of the current level of parents */
        while (parent < level_end) {
            IGRAPH_ALLOW_INTERRUPTION();
            for (j = 0; j < VECTOR(*branches)[k]; j++) {
                if (type == IGRAPH_TREE_IN) {
                    VECTOR(edges)[idx++] = child++;
                    VECTOR(edges)[idx++] = parent;
                } else {
                    VECTOR(edges)[idx++] = parent;
                    VECTOR(edges)[idx++] = child++;
                }
            }
            parent++;
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, type != IGRAPH_TREE_UNDIRECTED));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_regular_tree
 * \brief Creates a regular tree.
 *
 * All vertices of a regular tree, except its leaves, have the same total degree \p k.
 * This is different from a k-ary tree (\ref igraph_kary_tree()), where all
 * vertices have the same number of children, thus the degre of the root is
 * one less than the degree of the other internal vertices. Regular trees
 * are also referred to as Bethe lattices.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param h The height of the tree, i.e. the distance between the root and the leaves.
 * \param k The degree of the regular tree.
 * \param type Constant, gives whether to create a directed tree, and
 *        if this is the case, also its orientation. Possible values:
 *        \clist
 *        \cli IGRAPH_TREE_OUT
 *          directed tree, the edges point
 *          from the parents to their children.
 *        \cli IGRAPH_TREE_IN
 *          directed tree, the edges point from
 *          the children to their parents.
 *        \cli IGRAPH_TREE_UNDIRECTED
 *          undirected tree.
 *        \endclist
 *
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 *
 * \sa \ref igraph_kary_tree() to create k-ary tree where each vertex has the same
 * number of children, i.e. out-degree, instead of the same total degree.
 * \ref igraph_symmetric_tree() to use a different number of children at each level.
 *
 * \example examples/simple/igraph_regular_tree.c
 */

igraph_error_t igraph_regular_tree(igraph_t *graph, igraph_integer_t h, igraph_integer_t k, igraph_tree_mode_t type) {
    igraph_vector_int_t branching_counts;

    if (h < 1) {
        IGRAPH_ERRORF("Height of regular tree must be positive, got %" IGRAPH_PRId ".", IGRAPH_EINVAL, h);
    }
    if (k < 2 ) {
        IGRAPH_ERRORF("Degree of regular tree must be at least 2, got %" IGRAPH_PRId ".", IGRAPH_EINVAL, k);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&branching_counts, h);
    igraph_vector_int_fill(&branching_counts, k-1);
    if (h > 0) {
        VECTOR(branching_counts)[0] += 1;
    }

    IGRAPH_CHECK(igraph_symmetric_tree(graph, &branching_counts, type));

    igraph_vector_int_destroy(&branching_counts);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
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
igraph_error_t igraph_extended_chordal_ring(
    igraph_t *graph, igraph_integer_t nodes, const igraph_matrix_int_t *W,
    igraph_bool_t directed) {
    igraph_vector_int_t edges;
    igraph_integer_t period = igraph_matrix_int_ncol(W);
    igraph_integer_t nrow   = igraph_matrix_int_nrow(W);
    igraph_integer_t i, j, mpos = 0, epos = 0;

    if (nodes < 3) {
        IGRAPH_ERROR("An extended chordal ring has at least 3 nodes.", IGRAPH_EINVAL);
    }

    if (nodes % period != 0) {
        IGRAPH_ERROR("The period (number of columns in W) should divide the number of nodes.",
                     IGRAPH_EINVAL);
    }

    {
        /* ecount = nodes + nodes * nrow */
        igraph_integer_t no_of_edges2;
        IGRAPH_SAFE_MULT(nodes, nrow, &no_of_edges2);
        IGRAPH_SAFE_ADD(no_of_edges2, nodes, &no_of_edges2);
        IGRAPH_SAFE_MULT(no_of_edges2, 2, &no_of_edges2);
        IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges2);
    }

    for (i = 0; i < nodes - 1; i++) {
        VECTOR(edges)[epos++] = i;
        VECTOR(edges)[epos++] = i + 1;
    }
    VECTOR(edges)[epos++] = nodes - 1;
    VECTOR(edges)[epos++] = 0;

    if (nrow > 0) {
        for (i = 0; i < nodes; i++) {
            for (j = 0; j < nrow; j++) {
                igraph_integer_t offset = MATRIX(*W, j, mpos);
                igraph_integer_t v = (i + offset) % nodes;

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
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}
