
#include <igraph.h>

#include "bench.h"

int main(void) {
    igraph_t g;
    igraph_vector_int_list_t res;
    igraph_integer_t res_int;

    igraph_rng_seed(igraph_rng_default(), 42);
    BENCH_INIT();

    igraph_vector_int_list_init(&res, 0);

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 100, 3000, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    BENCH(" 1 Cliques in random graph with 100 vertices and 3000 edges",
          igraph_cliques(&g, &res, /* min_size= */ 0, /* max_size= */ 0);
         );
    igraph_destroy(&g);
    igraph_vector_int_list_clear(&res);

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 200, 10000, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    BENCH(" 2 Cliques in random graph with 200 vertices and 10000 edges, up to size 5",
          igraph_cliques(&g, &res, /* min_size= */ 0, /* max_size= */ 5);
         );
    igraph_vector_int_list_clear(&res);
    BENCH(" 3 Clique number of the same graph with 200 vertices and 10000 edges",
          igraph_clique_number(&g, &res_int);
         );
    igraph_vector_int_list_clear(&res);
    igraph_destroy(&g);

    igraph_vector_int_list_destroy(&res);

    return 0;
}
