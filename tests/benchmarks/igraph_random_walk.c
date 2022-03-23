#include <igraph.h>
#include "bench.h"

int main() {
    igraph_t graph;
    igraph_vector_int_t walk;
    igraph_vector_t weights;
    igraph_integer_t ec, i;

    igraph_rng_seed(igraph_rng_default(), 137);
    BENCH_INIT();

    igraph_vector_int_init(&walk, 0);
    igraph_vector_init(&weights, 0);

    /* create a small graph, and a compatible weight vector */
    igraph_de_bruijn(&graph, 3, 2); /* 9 vertices, 27 edges, average degree: 6 */
    ec = igraph_ecount(&graph);

    igraph_vector_resize(&weights, ec);
    for (i = 0; i < ec; ++i) {
        VECTOR(weights)[i] = igraph_rng_get_unif01(igraph_rng_default());
    }

    BENCH(" 1 Random edge walk,   directed,   unweighted, small graph  ",
          igraph_random_edge_walk(&graph, NULL, &walk, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 2 Random edge walk,   directed,   weighted,   small graph  ",
          igraph_random_edge_walk(&graph, &weights, &walk, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 3 Random vertex walk, directed,   unweighted, small graph  ",
          igraph_random_walk(&graph, &walk, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 4 Random vertex walk, directed,   weighted, small graph  ",
          igraph_random_walk(&graph, &weights, &walk, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 5 Random edge walk,   undirected, unweighted, small graph  ",
          igraph_random_edge_walk(&graph, NULL, &walk, 0, IGRAPH_ALL, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 6 Random edge walk,   undirected, weighted,   small graph  ",
          igraph_random_edge_walk(&graph, &weights, &walk, 0, IGRAPH_ALL, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 7 Random vertex walk, undirected, unweighted, small graph  ",
          igraph_random_walk(&graph, &walk, 0, IGRAPH_ALL, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 8 Random vertex walk, undirected, weighted, small graph  ",
          igraph_random_walk(&graph, &weights, &walk, 0, IGRAPH_ALL, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    igraph_destroy(&graph);

    /* create a big graph, and a compatible weight vector */
    igraph_de_bruijn(&graph, 8, 5); /* 32768 vertices, 262144 edges, average degree: 16 */
    ec = igraph_ecount(&graph);

    igraph_vector_resize(&weights, ec);
    for (i = 0; i < ec; ++i) {
        VECTOR(weights)[i] = igraph_rng_get_unif01(igraph_rng_default());
    }

    BENCH(" 9 Random edge walk,   directed,   unweighted, large graph  ",
          igraph_random_edge_walk(&graph, NULL, &walk, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 10 Random edge walk,   directed,   weighted,   large graph  ",
          igraph_random_edge_walk(&graph, &weights, &walk, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 11 Random vertex walk, directed,   unweighted, large graph  ",
          igraph_random_walk(&graph, &walk, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 12 Random vertex walk, directed,   weighted, large graph  ",
          igraph_random_walk(&graph, &weights, &walk, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );
    
    BENCH(" 13 Random edge walk,   undirected, unweighted, large graph  ",
          igraph_random_edge_walk(&graph, NULL, &walk, 0, IGRAPH_ALL, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 14 Random edge walk,   undirected, weighted,   large graph  ",
          igraph_random_edge_walk(&graph, &weights, &walk, 0, IGRAPH_ALL, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 15 Random vertex walk, undirected, unweighted, large graph  ",
          igraph_random_walk(&graph, &walk, 0, IGRAPH_ALL, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 16 Random vertex walk, undirected, weighted, large graph  ",
          igraph_random_walk(&graph, &weights, &walk, 0, IGRAPH_ALL, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 17 Short edge walk,    directed,   unweighted, large graph, x 100",
          REPEAT(igraph_random_edge_walk(&graph, NULL, &walk, 0, IGRAPH_OUT, 10000, IGRAPH_RANDOM_WALK_STUCK_RETURN), 100)
         );

    BENCH(" 18 Short edge walk,    directed,   weighted,   large graph, x 100",
          REPEAT(igraph_random_edge_walk(&graph, &weights, &walk, 0, IGRAPH_OUT, 10000, IGRAPH_RANDOM_WALK_STUCK_RETURN), 100)
         );

    BENCH(" 19 Short vertex walk,    directed,   unweighted, large graph, x 100",
          REPEAT(igraph_random_walk(&graph, NULL, &walk, 0, IGRAPH_OUT, 10000, IGRAPH_RANDOM_WALK_STUCK_RETURN), 100)
         );

    BENCH(" 20 Short vertex walk,    directed,   weighted,   large graph, x 100",
          REPEAT(igraph_random_walk(&graph, &weights, &walk, 0, IGRAPH_OUT, 10000, IGRAPH_RANDOM_WALK_STUCK_RETURN), 100)
         );

    BENCH(" 21 Short edge walk,    undirected, unweighted, large graph, x 100",
          REPEAT(igraph_random_edge_walk(&graph, NULL, &walk, 0, IGRAPH_ALL, 10000, IGRAPH_RANDOM_WALK_STUCK_RETURN), 100)
         );

    BENCH(" 22 Short edge walk,    undirected, weighted,   large graph, x 100",
          REPEAT(igraph_random_edge_walk(&graph, &weights, &walk, 0, IGRAPH_ALL, 10000, IGRAPH_RANDOM_WALK_STUCK_RETURN), 100)
         );
    
    BENCH(" 23 Short vertex walk,    undirected, unweighted, large graph, x 100",
          REPEAT(igraph_random_walk(&graph, NULL, &walk, 0, IGRAPH_ALL, 10000, IGRAPH_RANDOM_WALK_STUCK_RETURN), 100)
         );
    
    BENCH(" 24 Short vertex walk,    undirected, weighted,   large graph, x 100",
          REPEAT(igraph_random_walk(&graph, &weights, &walk, 0, IGRAPH_ALL, 10000, IGRAPH_RANDOM_WALK_STUCK_RETURN), 100)
         );

    igraph_destroy(&graph);

    igraph_vector_destroy(&weights);
    igraph_vector_int_destroy(&walk);

    return 0;
}
