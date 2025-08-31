/*
   igraph library.
   Copyright (C) 2024  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <igraph.h>
#include "bench.h"

int main(void) {
    igraph_t graph;
    igraph_vector_int_t vertices;
    igraph_vector_int_t edges;
    igraph_vector_t weights;
    igraph_int_t ec, i;

    igraph_rng_seed(igraph_rng_default(), 137);
    BENCH_INIT();

    igraph_vector_int_init(&vertices, 0);
    igraph_vector_int_init(&edges, 0);
    igraph_vector_init(&weights, 0);

    /* create a small graph, and a compatible weight vector */
    igraph_de_bruijn(&graph, 3, 2); /* 9 vertices, 27 edges, average degree: 6 */
    ec = igraph_ecount(&graph);

    igraph_vector_resize(&weights, ec);
    for (i = 0; i < ec; ++i) {
        VECTOR(weights)[i] = igraph_rng_get_unif01(igraph_rng_default());
    }

    /* Only edges */

    BENCH(" 1 Random edge walk,   directed,   unweighted, small graph  ",
          igraph_random_walk(&graph, NULL, NULL, &edges, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 2 Random edge walk,   directed,   weighted,   small graph  ",
          igraph_random_walk(&graph, &weights, NULL, &edges, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 3 Random edge walk,   undirected, unweighted, small graph  ",
        igraph_random_walk(&graph, NULL, NULL, &edges, 0, IGRAPH_ALL, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
    );

    BENCH(" 4 Random edge walk,   undirected, weighted,   small graph  ",
        igraph_random_walk(&graph, &weights, NULL, &edges, 0, IGRAPH_ALL, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
    );

    /* Only vertices */

    BENCH(" 5 Random vertex walk, directed,   unweighted, small graph  ",
          igraph_random_walk(&graph, NULL, &vertices, NULL, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 6 Random vertex walk, directed,   weighted, small graph  ",
          igraph_random_walk(&graph, &weights, &vertices, NULL, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 7 Random vertex walk, undirected, unweighted, small graph  ",
          igraph_random_walk(&graph, NULL, &vertices, NULL, 0, IGRAPH_ALL, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 8 Random vertex walk, undirected, weighted, small graph  ",
          igraph_random_walk(&graph, &weights, &vertices, NULL, 0, IGRAPH_ALL, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    /* Both edges and vertices */

    BENCH(" 9 Random walk,   directed,   unweighted, small graph  ",
        igraph_random_walk(&graph, NULL, &vertices, &edges, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
    );

    BENCH(" 10 Random walk,   directed,   weighted,   small graph  ",
        igraph_random_walk(&graph, &weights, &vertices, &edges, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
    );

    igraph_destroy(&graph);

    /* create a big graph, and a compatible weight vector */
    igraph_de_bruijn(&graph, 8, 5); /* 32768 vertices, 262144 edges, average degree: 16 */
    ec = igraph_ecount(&graph);

    igraph_vector_resize(&weights, ec);
    for (i = 0; i < ec; ++i) {
        VECTOR(weights)[i] = igraph_rng_get_unif01(igraph_rng_default());
    }

    /* Only edges */

    BENCH(" 11 Random edge walk,   directed,   unweighted, large graph  ",
          igraph_random_walk(&graph, NULL, NULL, &edges, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 12 Random edge walk,   directed,   weighted,   large graph  ",
          igraph_random_walk(&graph, &weights, NULL, &edges, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 13 Random edge walk,   undirected, unweighted, large graph  ",
        igraph_random_walk(&graph, NULL, NULL, &edges, 0, IGRAPH_ALL, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
    );

    BENCH(" 14 Random edge walk,   undirected, weighted,   large graph  ",
        igraph_random_walk(&graph, &weights, NULL, &edges, 0, IGRAPH_ALL, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
    );

    /* Only vertices */

    BENCH(" 15 Random vertex walk, directed,   unweighted, large graph  ",
          igraph_random_walk(&graph, NULL, &vertices, NULL, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 16 Random vertex walk, directed,   weighted, large graph  ",
          igraph_random_walk(&graph, &weights, &vertices, NULL, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 17 Random vertex walk, undirected, unweighted, large graph  ",
          igraph_random_walk(&graph, NULL, &vertices, NULL, 0, IGRAPH_ALL, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 18 Random vertex walk, undirected, weighted, large graph  ",
          igraph_random_walk(&graph, &weights, &vertices, NULL, 0, IGRAPH_ALL, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    /* Both edges and vertices */

    BENCH(" 19 Random walk,   directed,   unweighted, large graph  ",
        igraph_random_walk(&graph, NULL, &vertices, &edges, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
    );

    BENCH(" 20 Random walk,   directed,   weighted,   large graph  ",
        igraph_random_walk(&graph, &weights, &vertices, &edges, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
    );

    /* Only edges */

    BENCH(" 21 Short edge walk,    directed,   unweighted, large graph, x 100",
          REPEAT(igraph_random_walk(&graph, NULL, NULL, &edges, 0, IGRAPH_OUT, 10000, IGRAPH_RANDOM_WALK_STUCK_RETURN), 100)
         );

    BENCH(" 22 Short edge walk,    directed,   weighted,   large graph, x 100",
          REPEAT(igraph_random_walk(&graph, &weights, NULL, &edges, 0, IGRAPH_OUT, 10000, IGRAPH_RANDOM_WALK_STUCK_RETURN), 100)
         );

    BENCH(" 23 Short edge walk,    undirected, unweighted, large graph, x 100",
        REPEAT(igraph_random_walk(&graph, NULL, NULL, &edges, 0, IGRAPH_ALL, 10000, IGRAPH_RANDOM_WALK_STUCK_RETURN), 100)
    );

    BENCH(" 24 Short edge walk,    undirected, weighted,   large graph, x 100",
        REPEAT(igraph_random_walk(&graph, &weights, NULL, &edges, 0, IGRAPH_ALL, 10000, IGRAPH_RANDOM_WALK_STUCK_RETURN), 100)
    );

    /* Only vertices */

    BENCH(" 25 Short vertex walk,    directed,   unweighted, large graph, x 100",
          REPEAT(igraph_random_walk(&graph, NULL, &vertices, NULL, 0, IGRAPH_OUT, 10000, IGRAPH_RANDOM_WALK_STUCK_RETURN), 100)
         );

    BENCH(" 26 Short vertex walk,    directed,   weighted,   large graph, x 100",
          REPEAT(igraph_random_walk(&graph, &weights, &vertices, NULL, 0, IGRAPH_OUT, 10000, IGRAPH_RANDOM_WALK_STUCK_RETURN), 100)
         );

    BENCH(" 27 Short vertex walk,    undirected, unweighted, large graph, x 100",
          REPEAT(igraph_random_walk(&graph, NULL, &vertices, NULL, 0, IGRAPH_ALL, 10000, IGRAPH_RANDOM_WALK_STUCK_RETURN), 100)
         );

    BENCH(" 28 Short vertex walk,    undirected, weighted,   large graph, x 100",
          REPEAT(igraph_random_walk(&graph, &weights, &vertices, NULL, 0, IGRAPH_ALL, 10000, IGRAPH_RANDOM_WALK_STUCK_RETURN), 100)
         );

    /* Both edges and vertices */

    BENCH(" 29 Short random walk,    directed,   unweighted, large graph, x 100",
        REPEAT(igraph_random_walk(&graph, NULL, &vertices, &edges, 0, IGRAPH_OUT, 10000, IGRAPH_RANDOM_WALK_STUCK_RETURN), 100)
    );

    BENCH(" 30 Short random walk,    directed,   weighted,   large graph, x 100",
        REPEAT(igraph_random_walk(&graph, &weights, &vertices, &edges, 0, IGRAPH_OUT, 10000, IGRAPH_RANDOM_WALK_STUCK_RETURN), 100)
    );

    igraph_destroy(&graph);

    igraph_vector_destroy(&weights);
    igraph_vector_int_destroy(&vertices);
    igraph_vector_int_destroy(&edges);

    return 0;
}
