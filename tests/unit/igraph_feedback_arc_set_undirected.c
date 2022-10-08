
#include <igraph.h>

#include "test_utilities.h"

int main(void) {

    igraph_t graph;
    igraph_vector_int_t result;
    igraph_bool_t acyclic;

    igraph_rng_seed(igraph_rng_default(), 137);

    igraph_vector_int_init(&result, 0);

    /* Null graph */
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    igraph_feedback_arc_set(&graph, &result, NULL, IGRAPH_FAS_APPROX_EADES);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 0);
    igraph_destroy(&graph);

    /* Singleton graph */
    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    igraph_feedback_arc_set(&graph, &result, NULL, IGRAPH_FAS_APPROX_EADES);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 0);
    igraph_destroy(&graph);

    /* Two isolated vertices */
    igraph_empty(&graph, 2, IGRAPH_UNDIRECTED);
    igraph_feedback_arc_set(&graph, &result, NULL, IGRAPH_FAS_APPROX_EADES);
    IGRAPH_ASSERT(igraph_vector_int_size(&result) == 0);
    igraph_destroy(&graph);

    /* Single vertex with loop */
    acyclic = 0;
    igraph_small(&graph, 1, IGRAPH_UNDIRECTED,
                 0,0,
                 -1);
    igraph_feedback_arc_set(&graph, &result, NULL, IGRAPH_FAS_APPROX_EADES);
    igraph_delete_edges(&graph, igraph_ess_vector(&result));
    igraph_is_forest(&graph, &acyclic, NULL, IGRAPH_ALL);
    IGRAPH_ASSERT(acyclic);
    igraph_destroy(&graph);

    /* Random graph */
    acyclic = 0;
    igraph_erdos_renyi_game_gnm(&graph, 20, 40, IGRAPH_UNDIRECTED, IGRAPH_LOOPS);
    igraph_feedback_arc_set(&graph, &result, NULL, IGRAPH_FAS_APPROX_EADES);
    igraph_delete_edges(&graph, igraph_ess_vector(&result));
    igraph_is_forest(&graph, &acyclic, NULL, IGRAPH_ALL);
    IGRAPH_ASSERT(acyclic);
    igraph_destroy(&graph);

    igraph_vector_int_destroy(&result);

    VERIFY_FINALLY_STACK();

    return 0;
}
