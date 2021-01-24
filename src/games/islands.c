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
 * \param graph Pointer to an uninitialized graph object.
 * \param islands_n The number of islands in the graph.
 * \param islands_size The size of islands in the graph.
 * \param islands_pin The probability to create each possible edge into each island.
 * \param n_inter The number of edges to create between two islands.
 *
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid parameter
 *         \c IGRAPH_ENOMEM: there is not enough
 *         memory for the operation.
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 *
 */
igraph_long_t igraph_simple_interconnected_islands_game(
        igraph_t *graph,
        igraph_long_t islands_n,
        igraph_long_t islands_size,
        igraph_real_t islands_pin,
        igraph_long_t n_inter) {


    igraph_vector_t edges = IGRAPH_VECTOR_NULL;
    igraph_vector_t s = IGRAPH_VECTOR_NULL;
    igraph_long_t nbNodes;
    double maxpossibleedgesPerIsland;
    double maxedgesPerIsland;
    igraph_long_t nbEdgesInterIslands;
    double maxedges;
    igraph_long_t startIsland = 0;
    igraph_long_t endIsland = 0;
    igraph_long_t i, j, is;
    double myrand, last;
    igraph_long_t vsize;

    if (islands_n < 0) {
        IGRAPH_ERROR("Invalid number of islands", IGRAPH_EINVAL);
    }
    if (islands_size < 0) {
        IGRAPH_ERROR("Invalid size for islands", IGRAPH_EINVAL);
    }
    if (islands_pin < 0 || islands_pin > 1) {
        IGRAPH_ERROR("Invalid probability for islands", IGRAPH_EINVAL);
    }
    if ( (n_inter < 0) || (n_inter > islands_size) ) {
        IGRAPH_ERROR("Invalid number of inter-islands links", IGRAPH_EINVAL);
    }

    /* how much memory ? */
    nbNodes = islands_n * islands_size;
    maxpossibleedgesPerIsland = ((double)islands_size * ((double)islands_size - (double)1)) / (double)2;
    maxedgesPerIsland = islands_pin * maxpossibleedgesPerIsland;
    nbEdgesInterIslands = n_inter * (islands_n * (islands_n - 1)) / 2;
    maxedges = maxedgesPerIsland * islands_n + nbEdgesInterIslands;

    /* reserve enough space for all the edges */
    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&edges, (igraph_long_t) maxedges));

    RNG_BEGIN();

    /* first create all the islands */
    for (is = 1; is <= islands_n; is++) { /* for each island */

        /* index for start and end of nodes in this island */
        startIsland = islands_size * (is - 1);
        endIsland = startIsland + islands_size - 1;

        /* create the random numbers to be used (into s) */
        IGRAPH_VECTOR_INIT_FINALLY(&s, 0);
        IGRAPH_CHECK(igraph_vector_reserve(&s, (igraph_long_t) maxedgesPerIsland));

        last = RNG_GEOM(islands_pin);
        while (last < maxpossibleedgesPerIsland) { /* maxedgesPerIsland */
            IGRAPH_CHECK(igraph_vector_push_back(&s, last));
            myrand = RNG_GEOM(islands_pin);
            last += myrand; /* RNG_GEOM(islands_pin); */
            last += 1;
        }



        /* change this to edges ! */
        vsize = igraph_vector_size(&s);
        for (i = 0; i < vsize; i++) {
            igraph_long_t to = (igraph_long_t) floor((sqrt(8 * VECTOR(s)[i] + 1) + 1) / 2);
            igraph_long_t from = (igraph_long_t) (VECTOR(s)[i] - (((igraph_real_t)to) * (to - 1)) / 2);
            to += startIsland;
            from += startIsland;

            igraph_vector_push_back(&edges, from);
            igraph_vector_push_back(&edges, to);
        }

        /* clear the memory used for random number for this island */
        igraph_vector_destroy(&s);
        IGRAPH_FINALLY_CLEAN(1);


        /* create the links with other islands */
        for (i = is + 1; i <= islands_n; i++) { /* for each other island (not the previous ones) */

            for (j = 0; j < n_inter; j++) { /* for each link between islands */
                igraph_long_t from = (igraph_long_t) RNG_UNIF(startIsland, endIsland);
                igraph_long_t to = (igraph_long_t) RNG_UNIF((i - 1) * islands_size, i * islands_size);

                igraph_vector_push_back(&edges, from);
                igraph_vector_push_back(&edges, to);
            }

        }
    }

    RNG_END();

    /* actually fill the graph object */
    IGRAPH_CHECK(igraph_create(graph, &edges, nbNodes, 0));

    /* clean remaining things */
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
