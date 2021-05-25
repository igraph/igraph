
#include <igraph.h>

int main() {

    igraph_t g1, g2;
    igraph_bool_t iso;

    /* Seed the default random number generator and create a random graph. */

    igraph_rng_seed(igraph_rng_default(), 1122);

    igraph_erdos_renyi_game(&g1, IGRAPH_ERDOS_RENYI_GNP,
                            100, 3.0 / 100, /*directed=*/ 0, /*loops=*/ 0);

    /* Seed the generator with the same seed again,
     * and create a graph with the same method. */

    igraph_rng_seed(igraph_rng_default(), 1122);

    igraph_erdos_renyi_game(&g2, IGRAPH_ERDOS_RENYI_GNP,
                            100, 3.0 / 100, /*directed=*/ 0, /*loops=*/ 0);

    /* The two graphs will be identical. */

    igraph_is_same_graph(&g1, &g2, &iso);

    if (!iso) {
        return 1;
    }

    /* Destroy no longer needed data structures. */

    igraph_destroy(&g2);
    igraph_destroy(&g1);

    return 0;
}
