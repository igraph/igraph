#include <igraph.h>
#include <stdio.h>

int main(void) {
    igraph_integer_t nodes = 120, types = 4;

    igraph_matrix_t pref_matrix;
    igraph_matrix_init(&pref_matrix, types, types);

    igraph_rng_seed(igraph_rng_default(), 42);
    printf("Randomly generated graph with %" IGRAPH_PRId " nodes and %" IGRAPH_PRId " vertex types\n\n", nodes, types);

    /* Generate preference matrix giving connection probabilities for different vertex types */
    for (igraph_integer_t i = 0; i < types; i++) {
        for (igraph_integer_t j = 0; j < types; j++) {
            MATRIX(pref_matrix, i, j) = (i == j ? 0.1: 0.01);
        }
    }

    igraph_vector_int_t node_type_vec;
    igraph_vector_int_init(&node_type_vec, nodes);

    for (int i = 0; i < 5; i++) {
        igraph_real_t assortativity;
        igraph_t g;

        /* Generate undirected graph with 1000 nodes and 50 vertex types */
        igraph_preference_game(&g, nodes, types, /* type_dist= */ NULL, /* fixed_sizes= */ 1, &pref_matrix, &node_type_vec, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

        igraph_assortativity_nominal(&g, &node_type_vec, &assortativity, IGRAPH_UNDIRECTED, 1);
        printf("Assortativity before rewiring = %g\n", assortativity);

        /* Rewire graph */
        igraph_rewire(&g, 10 * igraph_ecount(&g), IGRAPH_REWIRING_SIMPLE);

        igraph_assortativity_nominal(&g, &node_type_vec, &assortativity, IGRAPH_UNDIRECTED, 1);
        printf("Assortativity after rewiring = %g\n\n", assortativity);

        igraph_destroy(&g);
    }
    igraph_vector_int_destroy(&node_type_vec);
    igraph_matrix_destroy(&pref_matrix);
}
