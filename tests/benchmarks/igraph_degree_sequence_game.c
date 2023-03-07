#include <igraph.h>

#include "bench.h"


int main(void) {
    igraph_t g, template;
    igraph_vector_int_t degrees;

    BENCH_INIT();
    igraph_rng_seed(igraph_rng_default(), 42);

    /* This benchmark indirectly tests the performance of igraph_set_t at the
     * moment because igraph_degree_sequence_game() (especially the
     * IGRAPH_DEGSEQ_CONFIGURATION_SIMPLE mode) heavily relies on sets.
     * This might change in future versions so the name of the benchmark
     * reflects the name of the function that is tested directly */

    igraph_vector_int_init(&degrees, 0);
    igraph_grg_game(&template, 1000, 0.04, /* torus = */ false, /* x = */ NULL, /* y = */ NULL);
    igraph_degree(&template, &degrees, igraph_vss_all(), IGRAPH_ALL, 1);
    BENCH(" 1 Degree sequence of GRG graph, N=1000, r=0.04, CONFIGURATION_SIMPLE",
          igraph_degree_sequence_game(&g, &degrees, /* indeg = */ 0, IGRAPH_DEGSEQ_CONFIGURATION_SIMPLE)
         );
    igraph_destroy(&g);
    igraph_destroy(&template);
    igraph_vector_int_destroy(&degrees);

    igraph_rng_seed(igraph_rng_default(), 42);
    igraph_vector_int_init(&degrees, 0);
    igraph_grg_game(&template, 10000, 0.013, /* torus = */ false, /* x = */ NULL, /* y = */ NULL);
    igraph_degree(&template, &degrees, igraph_vss_all(), IGRAPH_ALL, 1);
    BENCH(" 2 Degree sequence of GRG graph, N=10000, r=0.013, CONFIGURATION_SIMPLE",
          igraph_degree_sequence_game(&g, &degrees, /* indeg = */ 0, IGRAPH_DEGSEQ_CONFIGURATION_SIMPLE)
         );
    igraph_destroy(&g);
    igraph_destroy(&template);
    igraph_vector_int_destroy(&degrees);

    igraph_rng_seed(igraph_rng_default(), 42);
    igraph_vector_int_init(&degrees, 0);
    igraph_grg_game(&template, 10000, 0.0135, /* torus = */ false, /* x = */ NULL, /* y = */ NULL);
    igraph_degree(&template, &degrees, igraph_vss_all(), IGRAPH_ALL, 1);
    BENCH(" 3 Degree sequence of GRG graph, N=10000, r=0.0135, CONFIGURATION_SIMPLE",
          igraph_degree_sequence_game(&g, &degrees, /* indeg = */ 0, IGRAPH_DEGSEQ_CONFIGURATION_SIMPLE)
         );
    igraph_destroy(&g);
    igraph_destroy(&template);
    igraph_vector_int_destroy(&degrees);

    return 0;
}
