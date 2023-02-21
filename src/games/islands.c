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

#include "math/safe_intop.h"
#include "random/random_internal.h"

/**
 * \ingroup generators
 * \function igraph_simple_interconnected_islands_game
 * \brief Generates a random graph made of several interconnected islands, each island being a random graph.
 *
 * All islands are of the same size. Within an island, each edge is generated
 * with the same probability. A fixed number of additional edges are then
 * generated for each unordered pair of islands to connect them. The generated
 * graph is guaranteed to be simple.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param islands_n The number of islands in the graph.
 * \param islands_size The size of islands in the graph.
 * \param islands_pin The probability to create each possible edge within islands.
 * \param n_inter The number of edges to create between two islands. It may be
 *        larger than \p islands_size squared, but in this case it is assumed
 *        to be \p islands_size squared.
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
    igraph_integer_t start_index_of_island, start_index_of_other_island;
    igraph_integer_t i, j, is, from, to;
    igraph_real_t last;
    igraph_integer_t island_ecount;
    igraph_real_t nr_edges_reserved;

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

    number_of_inter_island_edges = islands_size * islands_size;
    if (n_inter > number_of_inter_island_edges) {
        IGRAPH_ERRORF(
            "Too many edges requested between islands, maximum possible "
            "is %" IGRAPH_PRId ", got %" IGRAPH_PRId ".",
            IGRAPH_EINVAL, number_of_inter_island_edges, n_inter
        );
    }

    /* how much memory ? */
    number_of_nodes = islands_n * islands_size;
    max_possible_edges_per_island = ((igraph_real_t)islands_size * ((igraph_real_t)islands_size - 1.0)) / 2.0;
    avg_edges_per_island = islands_pin * max_possible_edges_per_island;
    number_of_inter_island_edges = n_inter * (islands_n * (islands_n - 1)) / 2;

    nr_edges_reserved = 1.1 * avg_edges_per_island * islands_n + number_of_inter_island_edges;
    /* The cast of ECOUNT_MAX to double could change its value, which means in theory the size of
       the edges vector could still overflow, but only for very rare cases. */
    if (nr_edges_reserved > (double) (IGRAPH_ECOUNT_MAX ) || nr_edges_reserved > IGRAPH_MAX_EXACT_REAL) {
        IGRAPH_ERROR("Too many vertices, overflow in maximum number of edges.", IGRAPH_EOVERFLOW);
    }
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&edges, nr_edges_reserved * 2));

    IGRAPH_VECTOR_INIT_FINALLY(&s, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&s, 1.1 * avg_edges_per_island));

    RNG_BEGIN();

    /* first create all the islands */
    for (is = 0; is < islands_n; is++) { /* for each island */
        /* index for start and end of nodes in this island, both inclusive */
        start_index_of_island = islands_size * is;

        igraph_vector_clear(&s);

        last = RNG_GEOM(islands_pin);
        while (last < max_possible_edges_per_island) { /* avg_edges_per_island */
            IGRAPH_CHECK(igraph_vector_push_back(&s, last));
            last += RNG_GEOM(islands_pin);
            last += 1;
        }

        island_ecount = igraph_vector_size(&s);
        for (i = 0; i < island_ecount; i++) {
            to = floor((sqrt(8 * VECTOR(s)[i] + 1) + 1) / 2.0);
            from = VECTOR(s)[i] - (((igraph_real_t)to) * (to - 1)) / 2.0;
            to += start_index_of_island;
            from += start_index_of_island;

            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, from));
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to));
        }

        /* create the links with other islands */
        island_ecount = islands_size * islands_size;
        number_of_inter_island_edges = n_inter;
        for (i = is + 1; i < islands_n; i++) { /* for each other island (not the previous ones) */
            IGRAPH_CHECK(igraph_random_sample_real(&s, 0, island_ecount - 1, n_inter));

            start_index_of_other_island = i * islands_size;
            for (j = 0; j < n_inter; j++) { /* for each link between islands */
                from = VECTOR(s)[j] / islands_size;
                to = VECTOR(s)[j] - from * islands_size;
                from += start_index_of_island;
                to += start_index_of_other_island;

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
