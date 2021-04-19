
#include <igraph.h>

#include <time.h>

#include "test_utilities.inc"

static clock_t start;

/* Wait for at least a second before attempting interruption */
int interruption_handler(void *data) {
    if ( ((double) (clock() - start)) / CLOCKS_PER_SEC > 1.0 ) {
        IGRAPH_FINALLY_FREE();
        return IGRAPH_INTERRUPTED;
    } else {
        return IGRAPH_SUCCESS;
    }
}

int main() {
    igraph_t graph;
    igraph_vector_t res;
    igraph_error_handler_t *ehandler;

    igraph_vector_init(&res, 0);

    /* Skip test when igraph does not have GLPK support. */
    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1, -1);
    ehandler = igraph_set_error_handler(igraph_error_handler_ignore);
    if (igraph_feedback_arc_set(&graph, &res, NULL, IGRAPH_FAS_EXACT_IP) == IGRAPH_UNIMPLEMENTED) {
        igraph_destroy(&graph);
        igraph_vector_destroy(&res);
        return 77;
    }
    igraph_set_error_handler(ehandler);
    igraph_destroy(&graph);


    /* Current versions of GLPK will error if more than 100 million rows (MAX_M) are added.
       The graph size of 700 is chosen to just exceed this size. If future GLPK
       versions relax this restriction, the test will need to be updated accordinly. */
    igraph_full(&graph, 700, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);

    ehandler = igraph_set_error_handler(igraph_error_handler_printignore);
    IGRAPH_ASSERT(igraph_feedback_arc_set(&graph, &res, NULL, IGRAPH_FAS_EXACT_IP) == IGRAPH_EGLP);
    igraph_set_error_handler(ehandler);

    igraph_destroy(&graph);
    igraph_vector_destroy(&res);

    VERIFY_FINALLY_STACK();

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_vector_init(&res, 0);
    igraph_erdos_renyi_game(
                      &graph, IGRAPH_ERDOS_RENYI_GNM,
                      100, 200,
                      IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);

    igraph_set_interruption_handler(interruption_handler);
    ehandler = igraph_set_error_handler(igraph_error_handler_printignore);
    start = clock();
    IGRAPH_ASSERT(igraph_feedback_arc_set(&graph, &res, NULL, IGRAPH_FAS_EXACT_IP) == IGRAPH_INTERRUPTED);
    igraph_set_error_handler(ehandler);
    igraph_set_interruption_handler(NULL);

    igraph_destroy(&graph);

    igraph_vector_destroy(&res);

    VERIFY_FINALLY_STACK();

    return 0;
}
