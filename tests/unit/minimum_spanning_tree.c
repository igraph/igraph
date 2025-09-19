
#include <igraph.h>

#include "test_utilities.h"

void rand_weights(const igraph_t *graph, igraph_vector_t *weights) {
    igraph_int_t ecount = igraph_ecount(graph);
    igraph_vector_resize(weights, ecount);
    for (igraph_int_t i=0; i < ecount; i++) {
        VECTOR(*weights)[i] = RNG_UNIF01();
    }
}

igraph_real_t total_weight(const igraph_vector_int_t *edges, const igraph_vector_t *weights) {
    const igraph_int_t n = igraph_vector_int_size(edges);
    igraph_real_t sum = 0;
    for (igraph_int_t i=0; i < n; i++) {
        sum += VECTOR(*weights)[ VECTOR(*edges)[i] ];
    }
    return sum;
}

void check_forest(const igraph_t *g, const igraph_vector_int_t *edges) {
    igraph_t sg;
    igraph_bool_t is_forest;
    igraph_subgraph_from_edges(g, &sg, igraph_ess_vector(edges), false);
    igraph_is_forest(&sg, &is_forest, NULL, IGRAPH_ALL);
    IGRAPH_ASSERT(is_forest);
    igraph_destroy(&sg);
}

int main(void) {
    igraph_t g;
    igraph_vector_t weights;
    igraph_vector_int_t edges;
    igraph_int_t no_comps, no_nodes;

    igraph_rng_seed(igraph_rng_default(), 77685);

    igraph_vector_init(&weights, 0);
    igraph_vector_int_init(&edges, 0);

    igraph_erdos_renyi_game_gnm(&g, 50, 100, IGRAPH_UNDIRECTED, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED);
    rand_weights(&g, &weights);

    igraph_connected_components(&g, NULL, NULL, &no_comps, IGRAPH_WEAK);
    no_nodes = igraph_vcount(&g);

    printf("Automatic\n");
    igraph_minimum_spanning_tree(&g, &edges, &weights, IGRAPH_MST_AUTOMATIC);
    igraph_vector_int_sort(&edges);
    print_vector_int(&edges);
    printf("Total weight: %g\n", total_weight(&edges, &weights));
    check_forest(&g, &edges);
    IGRAPH_ASSERT(igraph_vector_int_size(&edges) == no_nodes - no_comps);

    printf("\nPrim\n");
    igraph_minimum_spanning_tree(&g, &edges, &weights, IGRAPH_MST_PRIM);
    igraph_vector_int_sort(&edges);
    print_vector_int(&edges);
    printf("Total weight: %g\n", total_weight(&edges, &weights));
    check_forest(&g, &edges);
    IGRAPH_ASSERT(igraph_vector_int_size(&edges) == no_nodes - no_comps);

    printf("\nKruskal\n");
    igraph_minimum_spanning_tree(&g, &edges, &weights, IGRAPH_MST_KRUSKAL);
    igraph_vector_int_sort(&edges);
    print_vector_int(&edges);
    printf("Total weight: %g\n", total_weight(&edges, &weights));
    check_forest(&g, &edges);
    IGRAPH_ASSERT(igraph_vector_int_size(&edges) == no_nodes - no_comps);

    printf("\nUnweighted\n");
    igraph_minimum_spanning_tree(&g, &edges, &weights, IGRAPH_MST_UNWEIGHTED);
    igraph_vector_int_sort(&edges);
    print_vector_int(&edges);
    check_forest(&g, &edges);
    IGRAPH_ASSERT(igraph_vector_int_size(&edges) == no_nodes - no_comps);

    igraph_destroy(&g);

    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&weights);

    VERIFY_FINALLY_STACK();

    return 0;
}
