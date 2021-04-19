
#include <igraph.h>

int main() {
    igraph_t graph;
    igraph_real_t result;

    /* Create a random preferential attachment graph. */
    igraph_barabasi_game(&graph, 30, /*power=*/ 1, 30, 0, 0, /*A=*/ 1,
                         IGRAPH_DIRECTED, IGRAPH_BARABASI_BAG,
                         /*start_from=*/ 0);

    /* Compute the average shortest path length. */
    igraph_average_path_length(&graph, &result, NULL, IGRAPH_UNDIRECTED, 1);
    printf("Average length of all-pairs shortest paths: %g\n", result);

    /* Destroy no-longer-needed objects. */
    igraph_destroy(&graph);

    return 0;
}
