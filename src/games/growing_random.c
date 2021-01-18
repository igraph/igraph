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
#include "igraph_random.h"

/**
 * \ingroup generators
 * \function igraph_growing_random_game
 * \brief Generates a growing random graph.
 *
 * </para><para>
 * This function simulates a growing random graph. We start out with
 * one vertex. In each step a new vertex is added and a number of new
 * edges are also added. These graphs are known to be different
 * from standard (not growing) random graphs.
 * \param graph Uninitialized graph object.
 * \param n The number of vertices in the graph.
 * \param m The number of edges to add in a time step (i.e. after
 *        adding a vertex).
 * \param directed Boolean, whether to generate a directed graph.
 * \param citation Boolean, if \c TRUE, the edges always
 *        originate from the most recently added vertex and are
 *        connected to a previous vertex.
 * \return Error code:
 *          \c IGRAPH_EINVAL: invalid
 *          \p n or \p m
 *          parameter.
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges.
 */
int igraph_growing_random_game(igraph_t *graph, igraph_integer_t n,
                               igraph_integer_t m, igraph_bool_t directed,
                               igraph_bool_t citation) {

    long int no_of_nodes = n;
    long int no_of_neighbors = m;
    long int no_of_edges;
    igraph_vector_t edges = IGRAPH_VECTOR_NULL;

    long int resp = 0;

    long int i, j;

    if (n < 0) {
        IGRAPH_ERROR("Invalid number of vertices", IGRAPH_EINVAL);
    }
    if (m < 0) {
        IGRAPH_ERROR("Invalid number of edges per step (m)", IGRAPH_EINVAL);
    }

    no_of_edges = (no_of_nodes - 1) * no_of_neighbors;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges * 2);

    RNG_BEGIN();

    for (i = 1; i < no_of_nodes; i++) {
        for (j = 0; j < no_of_neighbors; j++) {
            if (citation) {
                long int to = RNG_INTEGER(0, i - 1);
                VECTOR(edges)[resp++] = i;
                VECTOR(edges)[resp++] = to;
            } else {
                long int from = RNG_INTEGER(0, i);
                long int to = RNG_INTEGER(1, i);
                VECTOR(edges)[resp++] = from;
                VECTOR(edges)[resp++] = to;
            }
        }
    }

    RNG_END();

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}
