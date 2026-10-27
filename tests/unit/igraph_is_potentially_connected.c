#include <igraph.h>

#include "test_utilities.h"

static const char* bool_str(igraph_bool_t v) {
    if (v == true) {
        return " true";
    }
    else {
        return "false";
    }
}

/* Undirected case */
void potentially_connected_print_destroy(igraph_vector_int_t *ds) {
    igraph_bool_t g_simple, g_loops, g_multi, g_multiloops;
    igraph_bool_t c_simple, c_loops, c_multi, c_multiloops;

    print_vector_int(ds);

    igraph_is_graphical(ds, NULL, IGRAPH_MULTI_SW, &g_multi);
    igraph_is_graphical(ds, NULL, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW, &g_multiloops);
    igraph_is_graphical(ds, NULL, IGRAPH_LOOPS_SW, &g_loops);
    igraph_is_graphical(ds, NULL, IGRAPH_SIMPLE_SW, &g_simple);

    igraph_is_potentially_connected(ds, NULL, IGRAPH_SIMPLE_SW, IGRAPH_WEAK, &c_simple);
    igraph_is_potentially_connected(ds, NULL, IGRAPH_LOOPS_SW, IGRAPH_WEAK, &c_loops);
    igraph_is_potentially_connected(ds, NULL, IGRAPH_MULTI_SW, IGRAPH_WEAK, &c_multi);
    igraph_is_potentially_connected(ds, NULL, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW, IGRAPH_WEAK, &c_multiloops);

    printf("simple: %s, loops: %s, multi: %s, multiloops: %s\n\n",
           g_simple     ? bool_str(c_simple)     : "  n/a",
           g_loops      ? bool_str(c_loops)      : "  n/a",
           g_multi      ? bool_str(c_multi)      : "  n/a",
           g_multiloops ? bool_str(c_multiloops) : "  n/a");

    igraph_vector_int_destroy(ds);
}


/* Directed case */
void directed_potentially_connected_print_destroy(igraph_vector_int_t *ods, igraph_vector_int_t *ids) {
    igraph_bool_t g_simple, g_loops, g_multi, g_multiloops;
    igraph_bool_t c_simple, c_loops, c_multi, c_multiloops;

    print_vector_int(ods);
    print_vector_int(ids);

    igraph_is_graphical(ods, ids, IGRAPH_SIMPLE_SW, &g_simple);
    igraph_is_graphical(ods, ids, IGRAPH_LOOPS_SW, &g_loops);
    igraph_is_graphical(ods, ids, IGRAPH_MULTI_SW, &g_multi);
    igraph_is_graphical(ods, ids, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW, &g_multiloops);

    for (igraph_connectedness_t conn_mode = IGRAPH_WEAK; conn_mode <= IGRAPH_STRONG; conn_mode++) {

        igraph_is_potentially_connected(ods, ids, IGRAPH_SIMPLE_SW, conn_mode, &c_simple);
        if (conn_mode == IGRAPH_WEAK) {
            igraph_is_potentially_connected(ods, ids, IGRAPH_LOOPS_SW, conn_mode, &c_loops);
        } else {
            CHECK_ERROR(
                igraph_is_potentially_connected(ods, ids, IGRAPH_LOOPS_SW, conn_mode, &c_loops),
                IGRAPH_UNIMPLEMENTED
            );
        }
        igraph_is_potentially_connected(ods, ids, IGRAPH_MULTI_SW, conn_mode, &c_multi);
        igraph_is_potentially_connected(ods, ids, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW, conn_mode, &c_multiloops);

        printf("simple: %s, loops: %s, multi: %s, multiloops: %s\n",
               g_simple                                    ? bool_str(c_simple)     : "  n/a",
               (g_loops && conn_mode != IGRAPH_STRONG)     ? bool_str(c_loops)      : "  n/a",
               g_multi                                     ? bool_str(c_multi)      : "  n/a",
               g_multiloops                                ? bool_str(c_multiloops) : "  n/a");
    }
    printf("\n");

    igraph_vector_int_destroy(ods);
    igraph_vector_int_destroy(ids);
}


int main(void) {
    igraph_vector_int_t ds, ods, ids;

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
