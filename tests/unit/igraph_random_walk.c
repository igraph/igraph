
#include <igraph.h>

#include "test_utilities.inc"

int main() {
    igraph_t graph;
    igraph_vector_t walk, weights;
    igraph_integer_t ec, i;

    igraph_rng_seed(igraph_rng_default(), 137);

    igraph_vector_init(&walk, 0);
    igraph_vector_init(&weights, 0);

    /* This directed graph has loop edges.
       It also has multi-edges when considered as undirected. */
    igraph_de_bruijn(&graph, 3, 2);
    ec = igraph_ecount(&graph);

    /* unweighted, directed */
    igraph_random_edge_walk(&graph, NULL, &walk, 0, IGRAPH_OUT, 1000, IGRAPH_RANDOM_WALK_STUCK_RETURN);
    IGRAPH_ASSERT(igraph_vector_size(&walk) == 1000);

    /* unweighted, undirected */
    igraph_random_edge_walk(&graph, NULL, &walk, 0, IGRAPH_ALL, 1000, IGRAPH_RANDOM_WALK_STUCK_RETURN);
    IGRAPH_ASSERT(igraph_vector_size(&walk) == 1000);

    igraph_vector_resize(&weights, ec);
    for (i = 0; i < ec; ++i) {
        VECTOR(weights)[i] = igraph_rng_get_unif01(igraph_rng_default());
    }

    /* weighted, directed */
    igraph_random_edge_walk(&graph, &weights, &walk, 0, IGRAPH_OUT, 1000, IGRAPH_RANDOM_WALK_STUCK_RETURN);
    IGRAPH_ASSERT(igraph_vector_size(&walk) == 1000);

    /* weighted, undirecetd */
    igraph_random_edge_walk(&graph, &weights, &walk, 0, IGRAPH_ALL, 1000, IGRAPH_RANDOM_WALK_STUCK_RETURN);
    IGRAPH_ASSERT(igraph_vector_size(&walk) == 1000);

    igraph_destroy(&graph);

    /* 1-vertex graph, should get stuck */
    igraph_empty(&graph, 1, /* directed = */ 0);
    igraph_random_edge_walk(&graph, NULL, &walk, 0, IGRAPH_OUT, 1000, IGRAPH_RANDOM_WALK_STUCK_RETURN);
    IGRAPH_ASSERT(igraph_vector_size(&walk) == 0);
    igraph_destroy(&graph);

    igraph_vector_destroy(&weights);
    igraph_vector_destroy(&walk);

    VERIFY_FINALLY_STACK();

    return 0;
}
