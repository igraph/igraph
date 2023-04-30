#include <igraph.h>

#include "bench.h"
#include "../../src/core/indheap.h"
#include "igraph_types.h"

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define big 50000000

int main(void) {
    igraph_t graph;
    igraph_vector_t weights;

    BENCH_INIT();

    igraph_full(&graph, 1, IGRAPH_UNDIRECTED, 0);
    igraph_vector_init(&weights, igraph_ecount(&graph));
    BENCH("singleton graph", REPEAT(igraph_get_shortest_path_astar(&graph, NULL, NULL, 0, 0, &weights, IGRAPH_ALL, NULL, NULL), 1));
    igraph_vector_destroy(&weights);
    igraph_destroy(&graph);

    igraph_empty(&graph, big, IGRAPH_UNDIRECTED);
    igraph_vector_init(&weights, igraph_ecount(&graph));
    BENCH("graph with " TOSTRING(big) " nodes, no edges", REPEAT(igraph_get_shortest_path_astar(&graph, NULL, NULL, 0, 0, &weights, IGRAPH_ALL, NULL, NULL), 1));
    igraph_vector_destroy(&weights);
    igraph_destroy(&graph);

    igraph_ring(&graph, big, IGRAPH_UNDIRECTED, 0, 0);
    BENCH("vector init of " TOSTRING(big), igraph_vector_init(&weights, igraph_ecount(&graph)));
    BENCH("vector fill of " TOSTRING(big) " INFINITY", igraph_vector_fill(&weights, IGRAPH_INFINITY));
    BENCH("vector destroy", igraph_vector_destroy(&weights));
    igraph_vector_destroy(&weights);
    igraph_vector_init(&weights, igraph_ecount(&graph));
    BENCH("vector fill of " TOSTRING(big) " 1s", igraph_vector_fill(&weights, 1));

    BENCH("50000000 vertices ring graph astar from 0 to 0, 1x", REPEAT(igraph_get_shortest_path_astar(&graph, NULL, NULL, 0, 0, &weights, IGRAPH_ALL, NULL, NULL), 1));
    igraph_destroy(&graph);

    igraph_full(&graph, 10000, IGRAPH_UNDIRECTED, 0);
    igraph_vector_resize(&weights, igraph_ecount(&graph));
    igraph_vector_fill(&weights, 1);
    BENCH("10000 vertices full graph astar from 0 to 0, 1x", REPEAT(igraph_get_shortest_path_astar(&graph, NULL, NULL, 0, 0, &weights, IGRAPH_ALL, NULL, NULL), 1));

    igraph_vector_resize(&weights, big);
    igraph_vector_fill(&weights, 1);
    BENCH("vector_min on " TOSTRING(big) " numbers 1x", REPEAT(igraph_vector_min(&weights), 1));

    igraph_2wheap_t Q;
    BENCH("2wheap init of " TOSTRING(big) ", 1x", igraph_2wheap_init(&Q, big));

    igraph_integer_t *tmp;
    BENCH("IGRAPH_CALLOC of " TOSTRING(big) " igraph_integer_t, 1x", tmp = IGRAPH_CALLOC(big, igraph_integer_t));

    IGRAPH_FREE(tmp);
    BENCH("2wheap destroy", igraph_2wheap_destroy(&Q));
    igraph_vector_destroy(&weights);
    igraph_destroy(&graph);
}
