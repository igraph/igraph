
#include <igraph.h>

#include "bench.h"

void match(const igraph_t *graph,
           const igraph_t patt[], igraph_integer_t n,
           igraph_vector_int_list_t *maps) {

    igraph_vector_int_list_clear(maps);
    for (igraph_integer_t i=0; i < n; i++) {
        igraph_subisomorphic_lad(&patt[i], graph, NULL, NULL, NULL, maps, false, 0);
    }
}

#define NP 8

int main(void) {
    igraph_t graph;
    igraph_t patt[NP];
    igraph_vector_int_list_t maps;

    BENCH_INIT();

    igraph_vector_int_list_init(&maps, 0);

    for (igraph_integer_t i=0; i < NP; i++) {
        igraph_ring(&patt[i], i+1, IGRAPH_DIRECTED, false, true);
    }

    igraph_kautz(&graph, 3, 3);
    BENCH("1 Kautz(3,3) 10x", REPEAT(match(&graph, patt, NP, &maps), 10));
    igraph_destroy(&graph);

    igraph_kautz(&graph, 3, 4);
    BENCH("2 Kautz(3,4) 3x", REPEAT(match(&graph, patt, NP, &maps), 3));
    igraph_destroy(&graph);

    igraph_kautz(&graph, 4, 3);
    BENCH("3 Kautz(4,3) 1x", REPEAT(match(&graph, patt, NP, &maps), 1));
    igraph_destroy(&graph);

    for (igraph_integer_t i=0; i < NP; i++) {
        igraph_destroy(&patt[i]);
    }

    igraph_vector_int_list_destroy(&maps);

    return 0;
}
