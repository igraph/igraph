#include <igraph.h>

#include "bench.h"
#define TOSTR1(x) #x
#define TOSTR(x) TOSTR1(x)
#define BARBASI_GRAPH_VERTICE_COUNT 100000

void unit_test_subgraph_edges(
    const igraph_t *graph, igraph_vector_int_t *vertices,
    igraph_integer_t induced_subgraph_vertice_count, igraph_integer_t benchmark_count,
    const char* bench_message
) {
    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, 0);

    BENCH(bench_message,
          REPEAT(
              igraph_induced_subgraph_edges(graph, igraph_vss_range(0, induced_subgraph_vertice_count), &edges),
              benchmark_count
          )
         );

    igraph_vector_int_destroy(&edges);
}

void bench_induced_subgraph_edges(void) {
    igraph_t g;
    igraph_vector_int_t vertices;
    igraph_vector_int_t edges;
    igraph_vector_int_init_range(&vertices, 0, BARBASI_GRAPH_VERTICE_COUNT);
    igraph_vector_int_shuffle(&vertices);

    igraph_barabasi_game(&g, BARBASI_GRAPH_VERTICE_COUNT, 1, 100, NULL, true, 0,
                         IGRAPH_UNDIRECTED, IGRAPH_BARABASI_PSUMTREE_MULTIPLE, NULL);

#define REP 10

#define INDUCED_SUBGRAPH_VERTICE_COUNT 100

    unit_test_subgraph_edges(&g, &vertices, INDUCED_SUBGRAPH_VERTICE_COUNT, REP,
                             "benching induced_subgraph_edges for barabasi graph with vertices= \
" TOSTR(BARBASI_GRAPH_VERTICE_COUNT) " m=100, e` = " TOSTR(INDUCED_SUBGRAPH_VERTICE_COUNT) ",\
         " TOSTR(REP) " x"
                            );

#undef INDUCED_SUBGRAPH_VERTICE_COUNT

#define INDUCED_SUBGRAPH_VERTICE_COUNT 1000

    unit_test_subgraph_edges(&g, &vertices, INDUCED_SUBGRAPH_VERTICE_COUNT, REP,
                             "benching induced_subgraph_edges for barabasi graph with vertices= \
" TOSTR(BARBASI_GRAPH_VERTICE_COUNT) " m=100, e` = " TOSTR(INDUCED_SUBGRAPH_VERTICE_COUNT) ",\
         " TOSTR(REP) " x"
                            );

#undef INDUCED_SUBGRAPH_VERTICE_COUNT


#define INDUCED_SUBGRAPH_VERTICE_COUNT 10000

    unit_test_subgraph_edges(&g, &vertices, INDUCED_SUBGRAPH_VERTICE_COUNT, REP,
                             "benching induced_subgraph_edges for barabasi graph with vertices= \
" TOSTR(BARBASI_GRAPH_VERTICE_COUNT) " m=100, e` = " TOSTR(INDUCED_SUBGRAPH_VERTICE_COUNT) ",\
         " TOSTR(REP) " x"
                            );

#undef INDUCED_SUBGRAPH_VERTICE_COUNT

#undef REP


}

int main(void) {

    igraph_rng_seed(igraph_rng_default(), 137);
    BENCH_INIT();
    bench_induced_subgraph_edges();


    return 0;
}
