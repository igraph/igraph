#include <igraph.h>

#include "bench.h"

int main(void) {
    igraph_t graph;
    igraph_vector_t weights;

    BENCH_INIT();

    /*
    igraph_full(&graph, 10000, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_vector_init(&weights, igraph_ecount(&graph));
    igraph_vector_fill(&weights, 1);
    BENCH("singleton graph, with error checking", REPEAT(igraph_get_shortest_path_astar(&graph, NULL, NULL, 0, 0, &weights, IGRAPH_ALL, NULL, NULL), 100));

    igraph_vector_destroy(&weights);
    igraph_destroy(&graph);
    */

    igraph_ring(&graph, 50000000, IGRAPH_UNDIRECTED, 0, 0);
    igraph_vector_init(&weights, igraph_ecount(&graph));
    igraph_vector_fill(&weights, 1);
    BENCH("ring graph, with error checking", REPEAT(igraph_get_shortest_path_astar(&graph, NULL, NULL, 0, 0, &weights, IGRAPH_ALL, NULL, NULL), 1));

    igraph_vector_destroy(&weights);
    igraph_destroy(&graph);
}
