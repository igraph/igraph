
#include <igraph.h>

#include <time.h>

#include "test_utilities.h"

static clock_t start;

/* Wait for at least a second before attempting interruption */
igraph_bool_t interruption_handler(void) {
    if ( ((double) (clock() - start)) / CLOCKS_PER_SEC > 1.0 ) {
        IGRAPH_FINALLY_FREE();
        return true;
    } else {
        return false;
    }
}

int main(void) {
    igraph_t graph;
    igraph_vector_int_t res;
    igraph_error_handler_t *ehandler;

    igraph_vector_int_init(&res, 0);

    /* Skip test when igraph does not have GLPK support. */
    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1, -1);
    ehandler = igraph_set_error_handler(igraph_error_handler_ignore);
    if (igraph_feedback_arc_set(&graph, &res, NULL, IGRAPH_FAS_EXACT_IP_TI) == IGRAPH_UNIMPLEMENTED) {
        igraph_destroy(&graph);
        igraph_vector_int_destroy(&res);
        return 77;
    }
    igraph_set_error_handler(ehandler);
    igraph_destroy(&graph);


    /* Current versions of GLPK will error if more than 100 million rows (MAX_M) are added.
       The graph size of 700 is chosen to just exceed this size. If future GLPK
       versions relax this restriction, the test will need to be updated accordinly. */
    igraph_full(&graph, 700, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);

    ehandler = igraph_set_error_handler(igraph_error_handler_printignore);
    IGRAPH_ASSERT(igraph_feedback_arc_set(&graph, &res, NULL, IGRAPH_FAS_EXACT_IP_TI) == IGRAPH_FAILURE);
    igraph_set_error_handler(ehandler);

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&res);

    VERIFY_FINALLY_STACK();

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_vector_int_init(&res, 0);
    igraph_erdos_renyi_game_gnm(&graph, 100, 200, IGRAPH_DIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);

    igraph_set_interruption_handler(interruption_handler);
    ehandler = igraph_set_error_handler(igraph_error_handler_printignore);
    start = clock();
    IGRAPH_ASSERT(igraph_feedback_arc_set(&graph, &res, NULL, IGRAPH_FAS_EXACT_IP_TI) == IGRAPH_INTERRUPTED);
    igraph_set_error_handler(ehandler);
    igraph_set_interruption_handler(NULL);

    igraph_destroy(&graph);

    igraph_vector_int_destroy(&res);

    VERIFY_FINALLY_STACK();

    return 0;
}
