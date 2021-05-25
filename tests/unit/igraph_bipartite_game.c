
#include <igraph.h>

#include "test_utilities.inc"

int main() {
    igraph_t graph;
    igraph_bool_t bipartite;
    igraph_integer_t n1, n2, m;

    igraph_rng_seed(igraph_rng_default(), 947);

    n1 = 10; n2 = 20; m = 80;
    igraph_bipartite_game(&graph, NULL, IGRAPH_ERDOS_RENYI_GNM,
                          n1, n2, /* p */ 0, m,
                          IGRAPH_UNDIRECTED, IGRAPH_ALL);

    igraph_is_bipartite(&graph, &bipartite, NULL);

    IGRAPH_ASSERT(bipartite);
    IGRAPH_ASSERT(! igraph_is_directed(&graph));
    IGRAPH_ASSERT(igraph_vcount(&graph) == n1 + n2);
    IGRAPH_ASSERT(igraph_ecount(&graph) == m);

    igraph_destroy(&graph);

    n1 = 8; n2 = 15;
    igraph_bipartite_game(&graph, NULL, IGRAPH_ERDOS_RENYI_GNP,
                          n1, n2, 0.8, /* m */ 0,
                          IGRAPH_UNDIRECTED, IGRAPH_ALL);

    igraph_is_bipartite(&graph, &bipartite, NULL);

    IGRAPH_ASSERT(bipartite);
    IGRAPH_ASSERT(! igraph_is_directed(&graph));
    IGRAPH_ASSERT(igraph_vcount(&graph) == n1 + n2);
    IGRAPH_ASSERT(igraph_ecount(&graph) > 0); /* 0 is exceedingly unlikely */

    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
