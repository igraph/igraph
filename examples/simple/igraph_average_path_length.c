
#include <igraph.h>

int main(void) {
    igraph_t graph;
    igraph_real_t result;

    /* Initialize the library. */
    igraph_setup();

    /* Create a random preferential attachment graph. */
    igraph_barabasi_game(&graph, 30, /*power=*/ 1, 30, NULL, /*outpref=*/ false, /*A=*/ 1,
                         IGRAPH_DIRECTED, IGRAPH_BARABASI_BAG,
                         /*start_from=*/ NULL);

    /* Compute the average shortest path length. */
    igraph_average_path_length(&graph, /*weights=*/ NULL, &result,
                               /*unconn_pairs=*/ NULL, IGRAPH_UNDIRECTED, /*unconn=*/ true);
    printf("Average length of all-pairs shortest paths: %g\n", result);

    /* Destroy no-longer-needed objects. */
    igraph_destroy(&graph);

    return 0;
}
