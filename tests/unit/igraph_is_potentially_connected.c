#include <igraph.h>

#include "test_utilities.h"

/* Render a result value, flagging any value that was never written,
 * that is not 0 or 1, as UNSET so that a missing assignment is visible. */
static const char* bool_str(igraph_bool_t v) {
    if (v == 1) {
        return " true";
    }
    if (v == 0) {
        return "false";
    }
    return "UNSET";
}

/* Undirected case */
void potentially_connected_print_destroy(igraph_vector_int_t *ds) {
    int err;
    igraph_bool_t g_simple, g_loops, g_multi, g_multiloops;
    /* Initialized to a value that is neither 0 nor 1 so that a result the
     * function failed to write is reported as UNSET rather than as a
     * spurious true/false. */
    igraph_bool_t c_simple = 12345, c_loops = 12345, c_multi = 12345, c_multiloops = 12345;

    print_vector_int(ds);

    err = igraph_is_graphical(ds, NULL, IGRAPH_SIMPLE_SW, &g_simple);
    if (err != IGRAPH_SUCCESS) {
        printf("error!\n\n"); goto cleanup;
    }
    err = igraph_is_graphical(ds, NULL, IGRAPH_LOOPS_SW, &g_loops);
    if (err != IGRAPH_SUCCESS) {
        printf("error!\n\n"); goto cleanup;
    }
    err = igraph_is_graphical(ds, NULL, IGRAPH_MULTI_SW, &g_multi);
    if (err != IGRAPH_SUCCESS) {
        printf("error!\n\n"); goto cleanup;
    }
    err = igraph_is_graphical(ds, NULL, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW, &g_multiloops);
    if (err != IGRAPH_SUCCESS) {
        printf("error!\n\n"); goto cleanup;
    }
    
    if (g_simple) {
        err = igraph_is_potentially_connected(ds, NULL, IGRAPH_SIMPLE_SW, IGRAPH_WEAK, &c_simple);
        if (err != IGRAPH_SUCCESS) {
            printf("error!\n\n"); goto cleanup;
        }
    }
    if (g_loops) {
        err = igraph_is_potentially_connected(ds, NULL, IGRAPH_LOOPS_SW, IGRAPH_WEAK, &c_loops);
        if (err != IGRAPH_SUCCESS) {
            printf("error!\n\n"); goto cleanup;
        }
    }
    if (g_multi) {
        err = igraph_is_potentially_connected(ds, NULL, IGRAPH_MULTI_SW, IGRAPH_WEAK, &c_multi);
        if (err != IGRAPH_SUCCESS) {
            printf("error!\n\n"); goto cleanup;
        }
    }
    if (g_multiloops) {
        err = igraph_is_potentially_connected(ds, NULL, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW, IGRAPH_WEAK, &c_multiloops);
        if (err != IGRAPH_SUCCESS) {
            printf("error!\n\n"); goto cleanup;
        }
    }

    printf("simple: %s, loops: %s, multi: %s, multiloops: %s\n\n",
           g_simple     ? bool_str(c_simple)     : "  n/a",
           g_loops      ? bool_str(c_loops)      : "  n/a",
           g_multi      ? bool_str(c_multi)      : "  n/a",
           g_multiloops ? bool_str(c_multiloops) : "  n/a");
    fflush(stdout);
cleanup:
    igraph_vector_int_destroy(ds);
}


/* Directed case */
void directed_potentially_connected_print_destroy(igraph_vector_int_t *ods, igraph_vector_int_t *ids) {
    int err;
    igraph_bool_t g_simple, g_loops, g_multi, g_multiloops;
    /* Initialized to a value that is neither 0 nor 1 so that a result the
     * function failed to write is reported as UNSET rather than as a
     * spurious true/false. */
    igraph_bool_t c_simple = 12345, c_loops = 12345, c_multi = 12345, c_multiloops = 12345;

    print_vector_int(ods);
    print_vector_int(ids);

    err = igraph_is_graphical(ods, ids, IGRAPH_SIMPLE_SW, &g_simple);
    if (err != IGRAPH_SUCCESS) {
        printf("error!\n\n"); goto cleanup;
    }
    err = igraph_is_graphical(ods, ids, IGRAPH_LOOPS_SW, &g_loops);
    if (err != IGRAPH_SUCCESS) {
        printf("error!\n\n"); goto cleanup;
    }
    err = igraph_is_graphical(ods, ids, IGRAPH_MULTI_SW, &g_multi);
    if (err != IGRAPH_SUCCESS) {
        printf("error!\n\n"); goto cleanup;
    }
    err = igraph_is_graphical(ods, ids, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW, &g_multiloops);
    if (err != IGRAPH_SUCCESS) {
        printf("error!\n\n"); goto cleanup;
    }
    
    for (igraph_connectedness_t conn_mode = IGRAPH_WEAK; conn_mode <= IGRAPH_STRONG; conn_mode++) {
        if (g_simple) {
            err = igraph_is_potentially_connected(ods, ids, IGRAPH_SIMPLE_SW, conn_mode, &c_simple);
            if (err != IGRAPH_SUCCESS) {
                printf("error!\n\n"); goto cleanup;
            }
        }
        if (g_loops) {
            err = igraph_is_potentially_connected(ods, ids, IGRAPH_LOOPS_SW, conn_mode, &c_loops);
            if (err != IGRAPH_UNIMPLEMENTED && err != IGRAPH_SUCCESS) {
                printf("error!\n\n"); goto cleanup;
            }
        }
        if (g_multi) {
            err = igraph_is_potentially_connected(ods, ids, IGRAPH_MULTI_SW, conn_mode, &c_multi);
            if (err != IGRAPH_SUCCESS) {
                printf("error!\n\n"); goto cleanup;
            }
        }
        if (g_multiloops) {
            err = igraph_is_potentially_connected(ods, ids, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW, conn_mode, &c_multiloops);
            if (err != IGRAPH_SUCCESS) {
                printf("error!\n\n"); goto cleanup;
            }
        }
        
        printf("simple: %s, loops: %s, multi: %s, multiloops: %s\n",
               g_simple                                    ? bool_str(c_simple)     : "  n/a",
               (g_loops && conn_mode != IGRAPH_STRONG)     ? bool_str(c_loops)      : "  n/a",
               g_multi                                     ? bool_str(c_multi)      : "  n/a",
               g_multiloops                                ? bool_str(c_multiloops) : "  n/a");
    }
    printf("\n");
    fflush(stdout);

cleanup:
    igraph_vector_int_destroy(ods);
    igraph_vector_int_destroy(ids);
}


int main(void) {
    igraph_vector_int_t ds, ods, ids;

    igraph_set_error_handler(&igraph_error_handler_ignore);

    /* Undirected case: */

    /* Null graph */
    igraph_vector_int_init(&ds, 0);
    potentially_connected_print_destroy(&ds);
    
    /* Single vertex */
    igraph_vector_int_init_int_end(&ds, -1, 0, -1);
    potentially_connected_print_destroy(&ds);
    
    /* Single vertex with loop */
    igraph_vector_int_init_int_end(&ds, -1, 1, -1);
    potentially_connected_print_destroy(&ds);

    /* Single vertex many loops */
    igraph_vector_int_init_int_end(&ds, -1, 10, -1);
    potentially_connected_print_destroy(&ds);

    /* Two zeros */
    igraph_vector_int_init_int_end(&ds, -1, 0, 0, -1);
    potentially_connected_print_destroy(&ds);
    
    /* Short path */
    igraph_vector_int_init_int_end(&ds, -1, 1, 1, -1);
    potentially_connected_print_destroy(&ds);
        
    /* Should be disconnected for simple loopy */
    igraph_vector_int_init_int_end(&ds, -1, 2, 2, -1);
    potentially_connected_print_destroy(&ds);
        
    /* Should now be potentially connected for simple loopy */
    igraph_vector_int_init_int_end(&ds, -1, 4, 2, 2, -1);
    potentially_connected_print_destroy(&ds);

    /* Long path */
    igraph_vector_int_init_int_end(&ds, -1, 1, 2, 2, 2, 2, 2, 2, 2, 1, -1);
    potentially_connected_print_destroy(&ds);
    
    /* Cycle and isolated vertex */
    igraph_vector_int_init_int_end(&ds, -1, 2, 2, 2, 2, 2, 0, -1);
    potentially_connected_print_destroy(&ds);
    
    /* Tree */
    igraph_vector_int_init_int_end(&ds, -1, 3, 2, 2, 1, 1, 1, -1);
    potentially_connected_print_destroy(&ds);
    
    /* Almost tree */
    igraph_vector_int_init_int_end(&ds, -1, 3, 2, 1, 1, 1, 1, 1, -1);
    potentially_connected_print_destroy(&ds);

    /* Directed case: */

    /* Null graph */
    igraph_vector_int_init(&ods, 0);
    igraph_vector_int_init(&ids, 0);
    directed_potentially_connected_print_destroy(&ods, &ids);

    /* Single vertex */
    igraph_vector_int_init_int_end(&ods, -1, 0, -1);
    igraph_vector_int_init_int_end(&ids, -1, 0, -1);
    directed_potentially_connected_print_destroy(&ods, &ids);
    
    /* Single vertex with loop */
    igraph_vector_int_init_int_end(&ods, -1, 1, -1);
    igraph_vector_int_init_int_end(&ids, -1, 1, -1);
    directed_potentially_connected_print_destroy(&ods, &ids);
    
    /* Single vertex with many loops */
    igraph_vector_int_init_int_end(&ods, -1, 10, -1);
    igraph_vector_int_init_int_end(&ids, -1, 10, -1);
    directed_potentially_connected_print_destroy(&ods, &ids);

    /* Two isolated vertices */
    igraph_vector_int_init_int_end(&ods, -1, 0, 0, -1);
    igraph_vector_int_init_int_end(&ids, -1, 0, 0, -1);
    directed_potentially_connected_print_destroy(&ods, &ids);
    
    /* Single edge */
    igraph_vector_int_init_int_end(&ods, -1, 1, 0, -1);
    igraph_vector_int_init_int_end(&ids, -1, 0, 1, -1);
    directed_potentially_connected_print_destroy(&ods, &ids);
    
    /* 2-cycle */
    igraph_vector_int_init_int_end(&ods, -1, 1, 1, -1);
    igraph_vector_int_init_int_end(&ids, -1, 1, 1, -1);
    directed_potentially_connected_print_destroy(&ods, &ids);
    
    /* Exactly one strongly connected simple realization */
    igraph_vector_int_init_int_end(&ods, -1, 3, 2, 1, 1, -1);
    igraph_vector_int_init_int_end(&ids, -1, 1, 1, 3, 2, -1);
    directed_potentially_connected_print_destroy(&ods, &ids);
    
    /* Previous but add one edge to break strongly connected simple */
    igraph_vector_int_init_int_end(&ods, -1, 3, 3, 1, 1, -1);
    igraph_vector_int_init_int_end(&ids, -1, 1, 1, 3, 3, -1);
    directed_potentially_connected_print_destroy(&ods, &ids);


    VERIFY_FINALLY_STACK();

    return 0;
}
