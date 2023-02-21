
#include <igraph.h>

#include "bench.h"

int main(void) {
    igraph_t g;
    igraph_graph_list_t res;

    igraph_rng_seed(igraph_rng_default(), 42);
    BENCH_INIT();

    igraph_graph_list_init(&res, 0);

    igraph_empty(&g, 1000, IGRAPH_UNDIRECTED);
    BENCH(" 1 Decompose graph with 1000 isolated vertices",
          igraph_decompose(&g, &res, IGRAPH_WEAK, -1, -1);
         );
    igraph_destroy(&g);
    igraph_graph_list_clear(&res);

    igraph_empty(&g, 10000, IGRAPH_UNDIRECTED);
    BENCH(" 2 Decompose graph with 10000 isolated vertices",
          igraph_decompose(&g, &res, IGRAPH_WEAK, -1, -1);
         );
    igraph_destroy(&g);
    igraph_graph_list_clear(&res);

    igraph_empty(&g, 100000, IGRAPH_UNDIRECTED);
    BENCH(" 3 Decompose graph with 100000 isolated vertices",
          igraph_decompose(&g, &res, IGRAPH_WEAK, -1, -1);
         );
    igraph_destroy(&g);
    igraph_graph_list_clear(&res);

    igraph_graph_list_destroy(&res);

    return 0;
}
