
#include <igraph.h>

#include "bench.h"

void free_result(igraph_vector_ptr_t *res) {
    long int i, n;

    n = igraph_vector_ptr_size(res);
    for (i = 0; i < n; i++) {
        igraph_vector_t *v = VECTOR(*res)[i];
        igraph_vector_destroy(v);
        igraph_free(v);
    }

    igraph_vector_ptr_resize(res, 0);
}

int main() {
    igraph_t g;
    igraph_vector_ptr_t res;

    igraph_rng_seed(igraph_rng_default(), 42);
    BENCH_INIT();

    igraph_vector_ptr_init(&res, 0);

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 100, 3000, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    BENCH(" 1 Cliques in random graph with 100 vertices and 3000 edges",
          igraph_cliques(&g, &res, /* min_size= */ 0, /* max_size= */ 0);
         );
    igraph_destroy(&g);
    free_result(&res);

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 200, 10000, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    BENCH(" 2 Cliques in random graph with 200 vertices and 10000 edges, up to size 5",
          igraph_cliques(&g, &res, /* min_size= */ 0, /* max_size= */ 5);
         );
    igraph_destroy(&g);
    free_result(&res);

    igraph_vector_ptr_destroy(&res);

    return 0;
}
