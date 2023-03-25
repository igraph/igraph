#include <igraph.h>

#include "bench.h"


int main(void) {
    igraph_t g, template;
    igraph_vector_int_t degrees, outdeg, indeg;

    BENCH_INIT();
    igraph_rng_seed(igraph_rng_default(), 42);

    /* This benchmark indirectly tests the performance of igraph_set_t at the
     * moment because igraph_degree_sequence_game() (especially the
     * IGRAPH_DEGSEQ_CONFIGURATION_SIMPLE mode) heavily relies on sets.
     * This might change in future versions so the name of the benchmark
     * reflects the name of the function that is tested directly */

    igraph_vector_int_init(&degrees, 1000);
    igraph_vector_int_fill(&degrees, 7);
    BENCH(" 1 Degree sequence of undirected k-regular graph, N=1000, k=7, CONFIGURATION_SIMPLE",
          igraph_degree_sequence_game(&g, &degrees, /* indeg = */ NULL, IGRAPH_DEGSEQ_CONFIGURATION_SIMPLE)
         );
    igraph_destroy(&g);
    igraph_vector_int_destroy(&degrees);

    igraph_rng_seed(igraph_rng_default(), 42);
    igraph_vector_int_init(&degrees, 0);
    igraph_barabasi_game(&template, 1000, /* power = */ 1, /* m = */ 1,
        /* outseq = */ NULL, /* outpref = */ true, /* A = */ 0,
        IGRAPH_UNDIRECTED, IGRAPH_BARABASI_PSUMTREE, /* start_from = */ NULL
    );
    igraph_degree(&template, &degrees, igraph_vss_all(), IGRAPH_ALL, 1);
    BENCH(" 2 Degree sequence of undirected BA graph, N=1000, m=1, CONFIGURATION_SIMPLE",
          igraph_degree_sequence_game(&g, &degrees, /* indeg = */ NULL, IGRAPH_DEGSEQ_CONFIGURATION_SIMPLE)
         );
    igraph_destroy(&g);
    igraph_destroy(&template);
    igraph_vector_int_destroy(&degrees);

    igraph_rng_seed(igraph_rng_default(), 42);
    igraph_vector_int_init(&degrees, 0);
    igraph_erdos_renyi_game_gnm(&template, 200, 600, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_degree(&template, &degrees, igraph_vss_all(), IGRAPH_ALL, 1);
    BENCH(" 3 Degree sequence of undirected Erdos-Renyi graph, N=150, m=450, CONFIGURATION_SIMPLE",
          igraph_degree_sequence_game(&g, &degrees, /* indeg = */ NULL, IGRAPH_DEGSEQ_CONFIGURATION_SIMPLE)
         );
    igraph_destroy(&g);
    igraph_destroy(&template);
    igraph_vector_int_destroy(&degrees);

    igraph_rng_seed(igraph_rng_default(), 42);
    igraph_vector_int_init(&degrees, 0);
    igraph_grg_game(&template, 10000, 0.013, /* torus = */ false, /* x = */ NULL, /* y = */ NULL);
    igraph_degree(&template, &degrees, igraph_vss_all(), IGRAPH_ALL, true);
    BENCH(" 4 Degree sequence of GRG graph, N=10000, r=0.013, CONFIGURATION_SIMPLE",
          igraph_degree_sequence_game(&g, &degrees, /* indeg = */ NULL, IGRAPH_DEGSEQ_CONFIGURATION_SIMPLE)
         );
    igraph_destroy(&g);
    igraph_destroy(&template);
    igraph_vector_int_destroy(&degrees);

    igraph_vector_int_init(&degrees, 1000);
    igraph_vector_int_fill(&degrees, 5);
    BENCH(" 5 Degree sequence of directed k-regular graph, N=1000, k=5, CONFIGURATION_SIMPLE",
          igraph_degree_sequence_game(&g, &degrees, &degrees, IGRAPH_DEGSEQ_CONFIGURATION_SIMPLE)
         );
    igraph_destroy(&g);
    igraph_vector_int_destroy(&degrees);

    igraph_rng_seed(igraph_rng_default(), 42);
    igraph_vector_int_init(&outdeg, 0);
    igraph_vector_int_init(&indeg, 0);
    igraph_barabasi_game(&template, 5000, /* power = */ 1, /* m = */ 2,
        /* outseq = */ NULL, /* outpref = */ true, /* A = */ 0,
        IGRAPH_DIRECTED, IGRAPH_BARABASI_PSUMTREE, /* start_from = */ NULL
    );
    igraph_degree(&template, &outdeg, igraph_vss_all(), IGRAPH_OUT, true);
    igraph_degree(&template, &indeg, igraph_vss_all(), IGRAPH_IN, true);
    BENCH(" 6 Degree sequence of directed BA graph, N=500, m=2, CONFIGURATION_SIMPLE",
          igraph_degree_sequence_game(&g, &outdeg, &indeg, IGRAPH_DEGSEQ_CONFIGURATION_SIMPLE)
         );
    igraph_destroy(&g);
    igraph_destroy(&template);
    igraph_vector_int_destroy(&indeg);
    igraph_vector_int_destroy(&outdeg);

    igraph_rng_seed(igraph_rng_default(), 42);
    igraph_vector_int_init(&outdeg, 0);
    igraph_vector_int_init(&indeg, 0);
    igraph_erdos_renyi_game_gnm(&template, 15000, 45000, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    igraph_degree(&template, &outdeg, igraph_vss_all(), IGRAPH_OUT, true);
    igraph_degree(&template, &indeg, igraph_vss_all(), IGRAPH_IN, true);
    BENCH(" 7 Degree sequence of directed Erdos-Renyi graph, N=15000, m=45000, CONFIGURATION_SIMPLE",
          igraph_degree_sequence_game(&g, &outdeg, &indeg, IGRAPH_DEGSEQ_CONFIGURATION_SIMPLE)
         );
    igraph_destroy(&g);
    igraph_destroy(&template);
    igraph_vector_int_destroy(&indeg);
    igraph_vector_int_destroy(&outdeg);

    return 0;
}
