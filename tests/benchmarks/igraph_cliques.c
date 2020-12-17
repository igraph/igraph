
#include <igraph.h>

#include "bench.h"

int main() {
    igraph_t g;
    igraph_vector_ptr_t res;
    igraph_integer_t i, n;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 100, 3000, /* directed = */ 0, /* loops= */ 0);

    igraph_vector_ptr_init(&res, 0);

    BENCH("1 Cliques in random graph with 100 vertices and 3000 edges",
          igraph_cliques(&g, &res, /* min_size= */ 0, /* max_size= */ 0);
         );

    igraph_destroy(&g);

    n = igraph_vector_ptr_size(&res);
    for (i = 0; i < n; i++) {
        igraph_vector_t *v = VECTOR(res)[i];
        igraph_vector_destroy(v);
        igraph_free(v);
    }
    igraph_vector_ptr_destroy(&res);

    return 0;
}
