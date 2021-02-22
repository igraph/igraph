/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2003-2021 The igraph development team

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

#include "igraph_games.h"

#include "igraph_constructors.h"
#include "igraph_interface.h"
#include "igraph_random.h"

/* Uniform sampling of labelled trees (igraph_tree_game) */

/* The following implementation uniformly samples Prufer trees and converts
 * them to trees.
 */

static int igraph_i_tree_game_prufer(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed) {
    igraph_vector_int_t prufer;
    long i;

    if (directed) {
        IGRAPH_ERROR("The Prufer method for random tree generation does not support directed trees", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vector_int_init(&prufer, n - 2));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &prufer);

    RNG_BEGIN();

    for (i = 0; i < n - 2; ++i) {
        VECTOR(prufer)[i] = RNG_INTEGER(0, n - 1);
    }

    RNG_END();

    IGRAPH_CHECK(igraph_from_prufer(graph, &prufer));

    igraph_vector_int_destroy(&prufer);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/* The following implementation is based on loop-erased random walks and Wilson's algorithm
 * for uniformly sampling spanning trees. We effectively sample spanning trees of the complete
 * graph.
 */

/* swap two elements of a vector_int */
#define SWAP_INT_ELEM(vec, i, j) \
    { \
        igraph_integer_t temp; \
        temp = VECTOR(vec)[i]; \
        VECTOR(vec)[i] = VECTOR(vec)[j]; \
        VECTOR(vec)[j] = temp; \
    }

static int igraph_i_tree_game_loop_erased_random_walk(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed) {
    igraph_vector_t edges;
    igraph_vector_int_t vertices;
    igraph_vector_bool_t visited;
    long i, j, k;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 2 * (n - 1));

    IGRAPH_CHECK(igraph_vector_bool_init(&visited, n));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &visited);

    /* The vertices vector contains visited vertices between 0..k-1, unvisited ones between k..n-1. */
    IGRAPH_CHECK(igraph_vector_int_init_seq(&vertices, 0, n - 1));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &vertices);

    RNG_BEGIN();

    /* A simple implementation could be as below. This is for illustration only.
     * The actually implemented algorithm avoids unnecessary walking on the already visited
     * portion of the vertex set.
     */
    /*
    // pick starting point for the walk
    i = RNG_INTEGER(0, n-1);
    VECTOR(visited)[i] = 1;

    k=1;
    while (k < n) {
        // pick next vertex in the walk
        j = RNG_INTEGER(0, n-1);
        // if it has not been visited before, connect to the previous vertex in the sequence
        if (! VECTOR(visited)[j]) {
            VECTOR(edges)[2*k - 2] = i;
            VECTOR(edges)[2*k - 1] = j;
            VECTOR(visited)[j] = 1;
            k++;
        }
        i=j;
    }
    */

    i = RNG_INTEGER(0, n - 1);
    VECTOR(visited)[i] = 1;
    SWAP_INT_ELEM(vertices, 0, i);

    for (k = 1; k < n; ++k) {
        j = RNG_INTEGER(0, n - 1);
        if (VECTOR(visited)[VECTOR(vertices)[j]]) {
            i = VECTOR(vertices)[j];
            j = RNG_INTEGER(k, n - 1);
        }
        VECTOR(visited)[VECTOR(vertices)[j]] = 1;
        SWAP_INT_ELEM(vertices, k, j);
        VECTOR(edges)[2 * k - 2] = i;
        i = VECTOR(vertices)[k];
        VECTOR(edges)[2 * k - 1] = i;
    }

    RNG_END();

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));

    igraph_vector_int_destroy(&vertices);
    igraph_vector_bool_destroy(&visited);
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

#undef SWAP_INT_ELEM

/**
 * \ingroup generators
 * \function igraph_tree_game
 * \brief Generates a random tree with the given number of nodes.
 *
 * This function samples uniformly from the set of labelled trees,
 * i.e. it generates each labelled tree with the same probability.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param n The number of nodes in the tree.
 * \param directed Whether to create a directed tree. The edges are oriented away from the root.
 * \param method The algorithm to use to generate the tree. Possible values:
 *        \clist
 *        \cli IGRAPH_RANDOM_TREE_PRUFER
 *          This algorithm samples Pr&uuml;fer sequences uniformly, then converts them to trees.
 *          Directed trees are not currently supported.
 *        \cli IGRAPH_RANDOM_LERW
 *          This algorithm effectively performs a loop-erased random walk on the complete graph
 *          to uniformly sample its spanning trees (Wilson's algorithm).
 *        \endclist
 * \return Error code:
 *          \c IGRAPH_ENOMEM: there is not enough
 *           memory to perform the operation.
 *          \c IGRAPH_EINVAL: invalid tree size
 *
 * \sa \ref igraph_from_prufer()
 *
 */

int igraph_tree_game(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed, igraph_random_tree_t method) {
    if (n < 2) {
        IGRAPH_CHECK(igraph_empty(graph, n, directed));
        return IGRAPH_SUCCESS;
    }

    switch (method) {
    case IGRAPH_RANDOM_TREE_PRUFER:
        return igraph_i_tree_game_prufer(graph, n, directed);
    case IGRAPH_RANDOM_TREE_LERW:
        return igraph_i_tree_game_loop_erased_random_walk(graph, n, directed);
    default:
        IGRAPH_ERROR("Invalid method for random tree construction", IGRAPH_EINVAL);
    }
}
