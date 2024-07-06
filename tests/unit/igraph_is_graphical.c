
#include <igraph.h>

#include "test_utilities.h"

/* Undirected case */
void graphical_print_destroy(igraph_vector_int_t *ds) {
    int err;
    igraph_bool_t simple, loops, multi, multiloops;

    print_vector_int(ds);

    err = igraph_is_graphical(ds, NULL, IGRAPH_SIMPLE_SW, &simple);
    if (err != IGRAPH_SUCCESS) {
        printf("error!\n\n"); goto cleanup;
    }
    err = igraph_is_graphical(ds, NULL, IGRAPH_LOOPS_SW, &loops);
    if (err != IGRAPH_SUCCESS) {
        printf("error!\n\n"); goto cleanup;
    }
    err = igraph_is_graphical(ds, NULL, IGRAPH_MULTI_SW, &multi);
    if (err != IGRAPH_SUCCESS) {
        printf("error!\n\n"); goto cleanup;
    }
    err = igraph_is_graphical(ds, NULL, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW, &multiloops);
    if (err != IGRAPH_SUCCESS) {
        printf("error!\n\n"); goto cleanup;
    }

    printf("simple: %s, loops: %s, multi: %s, multiloops: %s\n\n",
           simple     ? " true" : "false",
           loops      ? " true" : "false",
           multi      ? " true" : "false",
           multiloops ? " true" : "false");

cleanup:
    igraph_vector_int_destroy(ds);
}


/* Directed case */
void digraphical_print_destroy(igraph_vector_int_t *ods, igraph_vector_int_t *ids) {
    int err;
    igraph_bool_t simple, loops, multi, multiloops;

    print_vector_int(ods);
    print_vector_int(ids);

    err = igraph_is_graphical(ods, ids, IGRAPH_SIMPLE_SW, &simple);
    if (err != IGRAPH_SUCCESS) {
        printf("error!\n\n"); goto cleanup;
    }
    err = igraph_is_graphical(ods, ids, IGRAPH_LOOPS_SW, &loops);
    if (err != IGRAPH_SUCCESS) {
        printf("error!\n\n"); goto cleanup;
    }
    err = igraph_is_graphical(ods, ids, IGRAPH_MULTI_SW, &multi);
    if (err != IGRAPH_SUCCESS) {
        printf("error!\n\n"); goto cleanup;
    }
    err = igraph_is_graphical(ods, ids, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW, &multiloops);
    if (err != IGRAPH_SUCCESS) {
        printf("error!\n\n"); goto cleanup;
    }

    printf("simple: %s, loops: %s, multi: %s, multiloops: %s\n\n",
           simple     ? " true" : "false",
           loops      ? " true" : "false",
           multi      ? " true" : "false",
           multiloops ? " true" : "false");

cleanup:
    igraph_vector_int_destroy(ods);
    igraph_vector_int_destroy(ids);
}


int main(void) {
    igraph_vector_int_t ds, ods, ids;

    igraph_set_error_handler(&igraph_error_handler_ignore);

    /* Undirected case: */

    /* Empty */
    igraph_vector_int_init(&ds, 0);
    graphical_print_destroy(&ds);

    /* All zeros */
    igraph_vector_int_init_int_end(&ds, -1, 0, 0, -1);
    graphical_print_destroy(&ds);

    igraph_vector_int_init_int_end(&ds, -1, 3, 3, 3, 3, 3, 3, 3, 3, -1);
    graphical_print_destroy(&ds);

    /* Undirected degree sequence with negative degree */
    igraph_vector_int_init_int_end(&ds, -1, 3, -2, 3, 3, 3, 3, 3, 3, -1);
    graphical_print_destroy(&ds);

    /* Undirected degree sequence with uneven sum */
    igraph_vector_int_init_int_end(&ds, -1, 3, 3, 3, 3, 3, 3, 3, -1);
    graphical_print_destroy(&ds);

    igraph_vector_int_init_int_end(&ds, -1, 4, 4, 5, 3, 6, 2, 2, 8, 1, 1, 10, -1);
    graphical_print_destroy(&ds);

    igraph_vector_int_init_int_end(&ds, -1, 3, 3, 2, 4, 1, 5, -1);
    graphical_print_destroy(&ds);

    igraph_vector_int_init_int_end(&ds, -1, 4, 7, 4, 7, 7, 8, 9, 9, 4, 6, 5, -1);
    graphical_print_destroy(&ds);

    igraph_vector_int_init_int_end(&ds, -1, 4, 4, 4, 4, 4, 1, 1, -1);
    graphical_print_destroy(&ds);

    igraph_vector_int_init_int_end(&ds, -1, 4, 4, 4, 4, 4, -1);
    graphical_print_destroy(&ds);

    igraph_vector_int_init_int_end(&ds, -1, 4, 4, 4, 4, 4, 4, 1, 1, -1);
    graphical_print_destroy(&ds);

    igraph_vector_int_init_int_end(&ds, -1, 3, 3, -1);
    graphical_print_destroy(&ds);

    igraph_vector_int_init_int_end(&ds, -1, 4, 4, 4, -1);
    graphical_print_destroy(&ds);

    igraph_vector_int_init_int_end(&ds, -1, 1, 2, 2, 3, -1);
    graphical_print_destroy(&ds);

    igraph_vector_int_init_int_end(&ds, -1, 1, 2, 3, -1);
    graphical_print_destroy(&ds);

    igraph_vector_int_init_int_end(&ds, -1, 1, 2, 5, -1);
    graphical_print_destroy(&ds);

    igraph_vector_int_init_int_end(&ds, -1, 1, 1, 4, -1);
    graphical_print_destroy(&ds);

    /* Extra cases for undirected simple with single loops */
    igraph_vector_int_init_int_end(&ds, -1, 7, 7, 3, 2, 2, 1, -1);
    graphical_print_destroy(&ds);

    igraph_vector_int_init_int_end(&ds, -1, 7, 7, 3, 3, 2, 2, -1);
    graphical_print_destroy(&ds);

    igraph_vector_int_init_int_end(&ds, -1, 6, 6, 6, 4, 2, 0, -1);
    graphical_print_destroy(&ds);

    igraph_vector_int_init_int_end(&ds, -1, 6, 6, 6, 4, 4, 0, -1);
    graphical_print_destroy(&ds);

    igraph_vector_int_init_int_end(&ds, -1, 6, 6, 6, 4, 3, 1, -1);
    graphical_print_destroy(&ds);

    /* The following two sequences are realizable as simple graphs.
     * The algorithm that checks this exits the last loop with these
     * two sequences. An earlier buggy version of the function failed
     * to set the result in this case. */
    igraph_vector_int_init_int_end(&ds, -1, 2, 2, 2, 2, 4, -1);
    graphical_print_destroy(&ds);

    igraph_vector_int_init_int_end(&ds, -1, 3, 0, 5, 3, 5, 3, 3, -1);
    graphical_print_destroy(&ds);


    /* Directed case: */

    /* Empty */
    igraph_vector_int_init(&ods, 0);
    igraph_vector_int_init(&ids, 0);
    digraphical_print_destroy(&ods, &ids);

    /* All zeros */
    igraph_vector_int_init_int_end(&ods, -1, 0, -1);
    igraph_vector_int_init_int_end(&ids, -1, 0, -1);
    digraphical_print_destroy(&ods, &ids);

    /* Different length, must throw an error */
    igraph_vector_int_init_int_end(&ods, -1, 1, 1, -1);
    igraph_vector_int_init_int_end(&ids, -1, 2, -1);
    digraphical_print_destroy(&ods, &ids);

    igraph_vector_int_init_int_end(&ods, -1, 0, 2, 3, 0, 4, 3, 1, 3, 4, 2, -1);
    igraph_vector_int_init_int_end(&ids, -1, 0, 3, 1, 3, 2, 4, 4, 1, 3, 1, -1);
    digraphical_print_destroy(&ods, &ids);

    igraph_vector_int_init_int_end(&ods, -1, 0, 2, 3, 0, 4, 3, 1, 3, 4, 2, -1);
    igraph_vector_int_init_int_end(&ids, -1, 0, 3, 1, -7, 2, 4, 4, 1, 3, 1, -1);
    digraphical_print_destroy(&ods, &ids);

    igraph_vector_int_init_int_end(&ods, -1, 0, 2, 3, 0, 4, 3, 1, 3, 4, 2, -1);
    igraph_vector_int_init_int_end(&ids, -1, 0, 3, 1, 2, 2, 4, 4, 1, 3, 1, -1);
    digraphical_print_destroy(&ods, &ids);

    igraph_vector_int_init_int_end(&ids, -1, 3, 3, 3, 3, 3, 3, 3, 3, 3, -1);
    igraph_vector_int_init_int_end(&ods, -1, 3, 3, 3, 3, 3, 3, 3, 3, 3, -1);
    digraphical_print_destroy(&ods, &ids);

    igraph_vector_int_init_int_end(&ids, -1, 1, 3, 2, 1, 3, 4, 3, 3, 1, 3, -1);
    igraph_vector_int_init_int_end(&ods, -1, 4, 1, 2, 3, 2, 3, 2, 3, 2, 2, -1);
    digraphical_print_destroy(&ods, &ids);

    igraph_vector_int_init_int_end(&ids, -1, 7, 4, 6, 4, 7, 8, 8, 8, 7, 4, -1);
    igraph_vector_int_init_int_end(&ods, -1, 8, 5, 6, 8, 6, 6, 5, 7, 5, 7, -1);
    digraphical_print_destroy(&ods, &ids);

    igraph_vector_int_init_int_end(&ids, -1, 3, 3, 1, 0, 2, 3, 0, 7, -1);
    igraph_vector_int_init_int_end(&ods, -1, 2, 2, 4, 3, 4, 3, 1, 0, -1);
    digraphical_print_destroy(&ods, &ids);

    /* Only one vertex with a non-zero out-degree. Regression test for bug #851 */
    igraph_vector_int_init_int_end(&ids, -1, 1, -1);
    igraph_vector_int_init_int_end(&ods, -1, 1, -1);
    digraphical_print_destroy(&ods, &ids);

    /* Another degree sequence when there is only
     * one vertex with a non-zero out-degree. Regression test for bug #851 */
    igraph_vector_int_init_int_end(&ids, -1, 2, 0, -1);
    igraph_vector_int_init_int_end(&ods, -1, 0, 2, -1);
    digraphical_print_destroy(&ods, &ids);

    igraph_vector_int_init_int_end(&ids, -1, 2, 2, -1);
    igraph_vector_int_init_int_end(&ods, -1, 2, 2, -1);
    digraphical_print_destroy(&ods, &ids);

    /* Valid directed graphical degree sequence. Regression test for bug #1092 */
    igraph_vector_int_init_int_end(&ids, -1, 1, 0, 1, -1);
    igraph_vector_int_init_int_end(&ods, -1, 0, 2, 0, -1);
    digraphical_print_destroy(&ods, &ids);

    /* Same as above, ids & ods exchanged. */
    igraph_vector_int_init_int_end(&ids, -1, 1, 0, 1, -1);
    igraph_vector_int_init_int_end(&ods, -1, 0, 2, 0, -1);
    digraphical_print_destroy(&ids, &ods);

    /* single loops: graphical, but multi-eges only: non-graphical */
    igraph_vector_int_init_int_end(&ids, -1, 1, 0, 2, -1);
    igraph_vector_int_init_int_end(&ods, -1, 0, 1, 2, -1);
    digraphical_print_destroy(&ids, &ods);

    /* Degree sequences of simple threshold digraphs:
     * These sequences make good test cases as they have a unique realization.
     * These sequnces were constructed based on point 3 of Theorem 1 in https://arxiv.org/abs/1212.1149
     * by computing a "closure" of a random digraphs.
     */
    igraph_vector_int_init_int_end(&ids, -1, 4, 3, 0, 3, 4, -1);
    igraph_vector_int_init_int_end(&ods, -1, 3, 3, 4, 3, 1, -1);
    digraphical_print_destroy(&ids, &ods);

    igraph_vector_int_init_int_end(&ids, -1, 4, 4, 3, 3, 4, -1);
    igraph_vector_int_init_int_end(&ods, -1, 4, 4, 4, 4, 2, -1);
    digraphical_print_destroy(&ids, &ods);

    igraph_vector_int_init_int_end(&ids, -1, 2, 4, 0, 3, 1, -1);
    igraph_vector_int_init_int_end(&ods, -1, 3, 2, 3, 1, 1, -1);
    digraphical_print_destroy(&ids, &ods);

    igraph_vector_int_init_int_end(&ids, -1, 11, 0, 0, 0, 7, 0, 9, 0, 13, 0, 0, 0, 0, 0, 0, -1);
    igraph_vector_int_init_int_end(&ods, -1, 3, 4, 4, 4, 3, 4, 3, 4, 2, 3, 2, 2, 1, 1, 0, -1);
    digraphical_print_destroy(&ids, &ods);

    igraph_vector_int_init_int_end(&ids, -1, 8, 4, 0, 0, 5, 0, 8, 9, 7, 0, 0, 11, 14, 13, 0, -1);
    igraph_vector_int_init_int_end(&ods, -1, 8, 8, 9, 9, 8, 8, 6, 5, 6, 4, 3, 2, 1, 1, 1, -1);
    digraphical_print_destroy(&ids, &ods);

    /* Some degree sequences with no simple realizations,
     * based on incrementing an in- and out-degree of the above. */
    igraph_vector_int_init_int_end(&ids, -1, 8, 4, 1, 0, 5, 0, 8, 9, 7, 0, 0, 11, 14, 13, 0, -1);
    igraph_vector_int_init_int_end(&ods, -1, 8, 8, 10, 9, 8, 8, 6, 5, 6, 4, 3, 2, 1, 1, 1, -1);
    digraphical_print_destroy(&ids, &ods);

    igraph_vector_int_init_int_end(&ids, -1, 9, 4, 0, 0, 5, 0, 8, 9, 7, 0, 0, 11, 14, 13, 0, -1);
    igraph_vector_int_init_int_end(&ods, -1, 8, 9, 9, 9, 8, 8, 6, 5, 6, 4, 3, 2, 1, 1, 1, -1);
    digraphical_print_destroy(&ids, &ods);

    igraph_vector_int_init_int_end(&ids, -1, 8, 4, 0, 0, 5, 0, 8, 9, 7, 0, 0, 11, 14, 14, 0, -1);
    igraph_vector_int_init_int_end(&ods, -1, 8, 8, 9, 9, 8, 8, 6, 5, 6, 4, 3, 3, 1, 1, 1, -1);
    digraphical_print_destroy(&ids, &ods);

    VERIFY_FINALLY_STACK();

    return 0;
}
