
#include <igraph.h>
#include <stdio.h>

#include "test_utilities.inc"

#define REALIZE1(ods, ids, et, method) \
    { \
        igraph_t graph; \
        int err; \
        err = igraph_realize_degree_sequence(&graph, ods, ids, et, method); \
        if (err == IGRAPH_SUCCESS) { \
            printf("\n"); \
            print_graph(&graph, stdout); \
            igraph_destroy(&graph); \
        } else if (err == IGRAPH_UNIMPLEMENTED) { \
            printf(" not implemented\n"); \
        } else { \
            printf(" not graphical\n"); \
        } \
    }

#define REALIZE2(ods, ids, et) \
    printf("Largest:"); \
    REALIZE1(ods, ids, et, IGRAPH_REALIZE_DEGSEQ_LARGEST); \
    printf("Smallest:"); \
    REALIZE1(ods, ids, et, IGRAPH_REALIZE_DEGSEQ_SMALLEST); \
    printf("Index:"); \
    REALIZE1(ods, ids, et, IGRAPH_REALIZE_DEGSEQ_INDEX);

#define UNDIRECTED_CHECK_DESTROY(ds) \
    print_vector_round(ds, stdout); \
    printf("\nSIMPLE GRAPH:\n"); \
    REALIZE2(ds, NULL, IGRAPH_SIMPLE_SW); \
    printf("\nLOOPLESS MULTIGRAPH:\n"); \
    REALIZE2(ds, NULL, IGRAPH_MULTI_SW); \
    printf("\nLOOPY MULTIGRAPH:\n"); \
    REALIZE2(ds, NULL, IGRAPH_MULTI_SW | IGRAPH_LOOPS_SW); \
    printf("\n\n"); \
    igraph_vector_destroy(ds);

#define DIRECTED_CHECK_DESTROY(ods, ids) \
    print_vector_round(ods, stdout); \
    print_vector_round(ids, stdout); \
    printf("\nSIMPLE GRAPH:\n"); \
    REALIZE2(ods, ids, IGRAPH_SIMPLE_SW); \
    printf("\nLOOPLESS MULTIGRAPH:\n"); \
    REALIZE2(ods, ids, IGRAPH_MULTI_SW); \
    printf("\nLOOPY MULTIGRAPH:\n"); \
    REALIZE2(ods, ids, IGRAPH_MULTI_SW | IGRAPH_LOOPS_SW); \
    printf("\n\n"); \
    igraph_vector_destroy(ods); \
    igraph_vector_destroy(ids);


int main() {
    igraph_vector_t ds, ods, ids;

    igraph_set_error_handler(&igraph_error_handler_ignore);

    /* Undirected */

    igraph_vector_init(&ds, 0);
    UNDIRECTED_CHECK_DESTROY(&ds);

    igraph_vector_init_int_end(&ds, -1, 1, 2, 2, 3, -1);
    UNDIRECTED_CHECK_DESTROY(&ds);

    /* contains negative degree */
    igraph_vector_init_int_end(&ds, -1, 1, 2, 2, -3, -1);
    UNDIRECTED_CHECK_DESTROY(&ds);

    /* odd sum */
    igraph_vector_init_int_end(&ds, -1, 1, 1, 2, 3, -1);
    UNDIRECTED_CHECK_DESTROY(&ds);

    igraph_vector_init_int_end(&ds, -1, 1, 2, 3, -1);
    UNDIRECTED_CHECK_DESTROY(&ds);

    igraph_vector_init_int_end(&ds, -1, 4, 4, 4, -1);
    UNDIRECTED_CHECK_DESTROY(&ds);

    igraph_vector_init_int_end(&ds, -1, 3, 5, -1);
    UNDIRECTED_CHECK_DESTROY(&ds);

    igraph_vector_init_int_end(&ds, -1, 5, 3, -1);
    UNDIRECTED_CHECK_DESTROY(&ds);

    igraph_vector_init_int_end(&ds, -1, 1, 3, 3, 4, 1, 2, 1, 1, 1, 3, -1);
    UNDIRECTED_CHECK_DESTROY(&ds);

    igraph_vector_init_int_end(&ds, -1, 2, 0, 3, 2, 2, 2, 2, 3, -1);
    UNDIRECTED_CHECK_DESTROY(&ds);

    /* Directed */

    igraph_vector_init(&ods, 0);
    igraph_vector_init(&ids, 0);
    DIRECTED_CHECK_DESTROY(&ods, &ids);

    igraph_vector_init_int_end(&ods, -1, 3, 0, 1, 1, 1, 1, 0, 1, -1);
    igraph_vector_init_int_end(&ids, -1, 2, 1, 0, 2, 2, 1, 0, 0, -1);
    DIRECTED_CHECK_DESTROY(&ods, &ids);

    igraph_vector_init_int_end(&ods, -1, 3, 1, 2, 3, 1, 2, 2, -1);
    igraph_vector_init_int_end(&ids, -1, 2, 2, 1, 2, 3, 2, 2, -1);
    DIRECTED_CHECK_DESTROY(&ods, &ids);


    VERIFY_FINALLY_STACK();

    return 0;
}
