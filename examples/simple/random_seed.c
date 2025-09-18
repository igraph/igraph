
#include <igraph.h>

int main(void) {
    igraph_t g1, g2;
    igraph_bool_t iso;

    /* Initialize the library. */
    igraph_setup();

    /* Seed the default random number generator and create a random graph. */
    igraph_rng_seed(igraph_rng_default(), 1122);

    igraph_erdos_renyi_game_gnp(&g1, 100, 3.0 / 100, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);

    /* Seed the generator with the same seed again,
     * and create a graph with the same method. */

    igraph_rng_seed(igraph_rng_default(), 1122);

    igraph_erdos_renyi_game_gnp(&g2, 100, 3.0 / 100, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);

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
