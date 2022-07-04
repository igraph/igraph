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
 * \function igraph_simple_interconnected_islands_game
 * \brief Generates a random graph made of several interconnected islands, each island being a random graph.
 *
 * All islands are of the same size. Within an island, each edge is generated
 * with the same probability. A fixed number of additional edges are then
 * generated for each unordered pair of islands to connect them. Note that
 * the islands themselves will \em not contain multiple edges, but there might
 * be multiple edges between the same pair of vertices between islands.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param islands_n The number of islands in the graph.
 * \param islands_size The size of islands in the graph.
 * \param islands_pin The probability to create each possible edge within islands.
 * \param n_inter The number of edges to create between two islands.
 *
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid parameter
 *         \c IGRAPH_ENOMEM: there is not enough memory for the operation.
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 *
 */
igraph_error_t igraph_simple_interconnected_islands_game(
        igraph_t *graph,
        igraph_integer_t islands_n,
        igraph_integer_t islands_size,
        igraph_real_t islands_pin,
        igraph_integer_t n_inter) {

    igraph_vector_int_t edges = IGRAPH_VECTOR_NULL;
    igraph_vector_t s = IGRAPH_VECTOR_NULL;
    igraph_integer_t number_of_nodes;
    igraph_real_t max_possible_edges_per_island;
    igraph_real_t avg_edges_per_island;
    igraph_integer_t number_of_inter_island_edges;
    igraph_integer_t start_island = 0;
    igraph_integer_t end_island = 0;
    igraph_integer_t i, j, is;
    igraph_real_t last;
    igraph_integer_t island_ecount;

    if (islands_n < 0) {
        IGRAPH_ERRORF("Number of islands cannot be negative, got %" IGRAPH_PRId ".", IGRAPH_EINVAL, islands_n);
    }
    if (islands_size < 0) {
        IGRAPH_ERRORF("Size of islands cannot be negative, got %" IGRAPH_PRId ".", IGRAPH_EINVAL, islands_size);
    }
    if (islands_pin < 0 || islands_pin > 1) {
        IGRAPH_ERRORF("Edge probability within islands should be between 0 and 1, got %g.", IGRAPH_EINVAL, islands_pin);
    }
    if (n_inter < 0) {
        IGRAPH_ERRORF("Number of inter-island links cannot be negative, got %" IGRAPH_PRId ".", IGRAPH_EINVAL, n_inter);
    }

    /* how much memory ? */
    number_of_nodes = islands_n * islands_size;
    max_possible_edges_per_island = ((igraph_real_t)islands_size * ((igraph_real_t)islands_size - 1.0)) / 2.0;
    avg_edges_per_island = islands_pin * max_possible_edges_per_island;
    number_of_inter_island_edges = n_inter * (islands_n * (islands_n - 1)) / 2;

    /* reserve enough space for all the edges */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&edges, 1.1 * avg_edges_per_island * islands_n + number_of_inter_island_edges));

    IGRAPH_VECTOR_INIT_FINALLY(&s, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&s, 1.1 * avg_edges_per_island));

    RNG_BEGIN();

    /* first create all the islands */
    for (is = 0; is < islands_n; is++) { /* for each island */

        /* index for start and end of nodes in this island */
        start_island = islands_size * is;
        end_island = start_island + islands_size - 1;

        igraph_vector_clear(&s);

        last = RNG_GEOM(islands_pin);
        while (last < max_possible_edges_per_island) { /* avg_edges_per_island */
            IGRAPH_CHECK(igraph_vector_push_back(&s, last));
            last += RNG_GEOM(islands_pin);
            last += 1;
        }

        island_ecount = igraph_vector_size(&s);
        for (i = 0; i < island_ecount; i++) {
            igraph_integer_t to = floor((sqrt(8 * VECTOR(s)[i] + 1) + 1) / 2.0);
            igraph_integer_t from = VECTOR(s)[i] - (((igraph_real_t)to) * (to - 1)) / 2.0;
            to += start_island;
            from += start_island;

            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, from));
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to));
        }

        /* create the links with other islands */
        for (i = is + 1; i < islands_n; i++) { /* for each other island (not the previous ones) */

            for (j = 0; j < n_inter; j++) { /* for each link between islands */
                igraph_integer_t from = RNG_INTEGER(start_island, end_island);
                igraph_integer_t to = RNG_INTEGER(i * islands_size, (i + 1) * islands_size - 1);

                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, from));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to));
            }

        }
    }

    igraph_vector_destroy(&s);
    IGRAPH_FINALLY_CLEAN(1);

    RNG_END();

    /* actually fill the graph object */
    IGRAPH_CHECK(igraph_create(graph, &edges, number_of_nodes, 0));

    /* clean remaining things */
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
