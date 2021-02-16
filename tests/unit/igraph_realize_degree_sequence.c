
#include <igraph.h>
#include <stdio.h>

#include "test_utilities.inc"

void realize1(igraph_vector_t *ods, igraph_vector_t *ids, igraph_edge_type_sw_t et, igraph_realize_degseq_t method) {
    igraph_t graph;
    int err;
    err = igraph_realize_degree_sequence(&graph, ods, ids, et, method);
    if (err == IGRAPH_SUCCESS) {
        printf("\n");
        print_graph(&graph);
        igraph_destroy(&graph);
    } else if (err == IGRAPH_UNIMPLEMENTED) {
        printf(" not implemented\n");
    } else {
        printf(" not graphical\n");
    }
}

void realize2(igraph_vector_t *ods, igraph_vector_t *ids, igraph_edge_type_sw_t et) {
    printf("Largest:");
    realize1(ods, ids, et, IGRAPH_REALIZE_DEGSEQ_LARGEST);
    printf("Smallest:");
    realize1(ods, ids, et, IGRAPH_REALIZE_DEGSEQ_SMALLEST);
    printf("Index:");
    realize1(ods, ids, et, IGRAPH_REALIZE_DEGSEQ_INDEX);
}

void undirected_print_destroy(igraph_vector_t *ds) {
    print_vector_round(ds);
    printf("\nSIMPLE GRAPH:\n");
    realize2(ds, NULL, IGRAPH_SIMPLE_SW);
    printf("\nLOOPLESS MULTIGRAPH:\n");
    realize2(ds, NULL, IGRAPH_MULTI_SW);
    printf("\nLOOPY MULTIGRAPH:\n");
    realize2(ds, NULL, IGRAPH_MULTI_SW | IGRAPH_LOOPS_SW);
    printf("\n\n");
    igraph_vector_destroy(ds);
}

void directed_print_destroy(igraph_vector_t *ods, igraph_vector_t *ids) {
    print_vector_round(ods);
    print_vector_round(ids);
    printf("\nSIMPLE GRAPH:\n");
    realize2(ods, ids, IGRAPH_SIMPLE_SW);
    printf("\nLOOPLESS MULTIGRAPH:\n");
    realize2(ods, ids, IGRAPH_MULTI_SW);
    printf("\nLOOPY MULTIGRAPH:\n");
    realize2(ods, ids, IGRAPH_MULTI_SW | IGRAPH_LOOPS_SW);
    printf("\n\n");
    igraph_vector_destroy(ods);
    igraph_vector_destroy(ids);
}


int main() {
    igraph_vector_t ds, ods, ids;

    igraph_set_error_handler(&igraph_error_handler_ignore);

    /* Undirected */

    igraph_vector_init(&ds, 0);
    undirected_print_destroy(&ds);

    igraph_vector_init_int_end(&ds, -1, 1, 2, 2, 3, -1);
    undirected_print_destroy(&ds);

    /* contains negative degree */
    igraph_vector_init_int_end(&ds, -1, 1, 2, 2, -3, -1);
    undirected_print_destroy(&ds);

    /* odd sum */
    igraph_vector_init_int_end(&ds, -1, 1, 1, 2, 3, -1);
    undirected_print_destroy(&ds);

    igraph_vector_init_int_end(&ds, -1, 1, 2, 3, -1);
    undirected_print_destroy(&ds);

    igraph_vector_init_int_end(&ds, -1, 4, 4, 4, -1);
    undirected_print_destroy(&ds);

    igraph_vector_init_int_end(&ds, -1, 3, 5, -1);
    undirected_print_destroy(&ds);

    igraph_vector_init_int_end(&ds, -1, 5, 3, -1);
    undirected_print_destroy(&ds);

    igraph_vector_init_int_end(&ds, -1, 1, 3, 3, 4, 1, 2, 1, 1, 1, 3, -1);
    undirected_print_destroy(&ds);

    igraph_vector_init_int_end(&ds, -1, 2, 0, 3, 2, 2, 2, 2, 3, -1);
    undirected_print_destroy(&ds);

    /* Directed */

    igraph_vector_init(&ods, 0);
    igraph_vector_init(&ids, 0);
    directed_print_destroy(&ods, &ids);

    igraph_vector_init_int_end(&ods, -1, 3, 0, 1, 1, 1, 1, 0, 1, -1);
    igraph_vector_init_int_end(&ids, -1, 2, 1, 0, 2, 2, 1, 0, 0, -1);
    directed_print_destroy(&ods, &ids);

    igraph_vector_init_int_end(&ods, -1, 3, 1, 2, 3, 1, 2, 2, -1);
    igraph_vector_init_int_end(&ids, -1, 2, 2, 1, 2, 3, 2, 2, -1);
    directed_print_destroy(&ods, &ids);

    /* single loops: graphical, but multi-eges only: non-graphical */
    igraph_vector_init_int_end(&ids, -1, 1, 0, 2, -1);
    igraph_vector_init_int_end(&ods, -1, 0, 1, 2, -1);
    directed_print_destroy(&ods, &ids);

    /* same as before, different ordering, to test the "index" method */
    igraph_vector_init_int_end(&ids, -1, 2, 0, 1, -1);
    igraph_vector_init_int_end(&ods, -1, 2, 1, 0, -1);
    directed_print_destroy(&ods, &ids);

    /* same as before, different ordering, to test the "index" method */
    igraph_vector_init_int_end(&ids, -1, 0, 2, 1, -1);
    igraph_vector_init_int_end(&ods, -1, 1, 2, 0, -1);
    directed_print_destroy(&ods, &ids);

    igraph_vector_init_int_end(&ids, -1, 2, 0, -1);
    igraph_vector_init_int_end(&ods, -1, 0, 2, -1);
    directed_print_destroy(&ods, &ids);

    /* simple complete graph on 4 vertices */
    igraph_vector_init_int_end(&ids, -1, 3, 3, 3, 3, -1);
    igraph_vector_init_int_end(&ods, -1, 3, 3, 3, 3, -1);
    directed_print_destroy(&ods, &ids);

    VERIFY_FINALLY_STACK();

    return 0;
}
