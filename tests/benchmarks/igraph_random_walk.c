
#include <igraph.h>
#include "bench.h"

int main() {
    igraph_t graph;
    igraph_vector_t walk, weights;
    igraph_integer_t ec, i;

    igraph_rng_seed(igraph_rng_default(), 137);
    BENCH_INIT();

    igraph_vector_init(&walk, 0);
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

    igraph_to_undirected(&graph, IGRAPH_TO_UNDIRECTED_EACH, NULL);

    BENCH(" 4 Random edge walk,   undirected, unweighted, small graph  ",
          igraph_random_edge_walk(&graph, NULL, &walk, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 5 Random edge walk,   undirected, weighted,   small graph  ",
          igraph_random_edge_walk(&graph, &weights, &walk, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 6 Random vertex walk, undirected, unweighted, small graph  ",
          igraph_random_walk(&graph, &walk, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    igraph_destroy(&graph);

    /* create a big graph, and a compatible weight vector */
    igraph_de_bruijn(&graph, 8, 5); /* 32768 vertices, 262144 edges, average degree: 16 */
    ec = igraph_ecount(&graph);

    igraph_vector_resize(&weights, ec);
    for (i = 0; i < ec; ++i) {
        VECTOR(weights)[i] = igraph_rng_get_unif01(igraph_rng_default());
    }

    BENCH(" 7 Random edge walk,   directed,   unweighted, large graph  ",
          igraph_random_edge_walk(&graph, NULL, &walk, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 8 Random edge walk,   directed,   weighted,   large graph  ",
          igraph_random_edge_walk(&graph, &weights, &walk, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH(" 9 Random vertex walk, directed,   unweighted, large graph  ",
          igraph_random_walk(&graph, &walk, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    igraph_to_undirected(&graph, IGRAPH_TO_UNDIRECTED_EACH, NULL);

    BENCH("10 Random edge walk,   undirected, unweighted, large graph  ",
          igraph_random_edge_walk(&graph, NULL, &walk, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH("11 Random edge walk,   undirected, weighted,   large graph  ",
          igraph_random_edge_walk(&graph, &weights, &walk, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    BENCH("12 Random vertex walk, undirected, unweighted, large graph  ",
          igraph_random_walk(&graph, &walk, 0, IGRAPH_OUT, 50000000, IGRAPH_RANDOM_WALK_STUCK_RETURN)
         );

    igraph_destroy(&graph);

    igraph_vector_destroy(&weights);
    igraph_vector_destroy(&walk);

    return 0;
}
