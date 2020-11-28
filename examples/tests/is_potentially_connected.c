
#include <igraph.h>

#include "test_utilities.inc"

void pc_undirected_print_destroy(igraph_vector_t *ds) {
    igraph_bool_t res;
    int err;

    printf("\n");
    print_vector_round(ds, stdout);
    err = igraph_is_potentially_connected_degree_sequence(ds, NULL, &res);
    if (err != IGRAPH_SUCCESS) {
        printf("error!\n"); goto cleanup;
    }
    printf("%s\n", res ? "true" : "false");

cleanup:
    igraph_vector_destroy(ds);
}

void pc_directed_print_destroy(igraph_vector_t *ods, igraph_vector_t *ids) {
    igraph_bool_t res;
    int err;

    printf("\n");
    print_vector_round(ods, stdout);
    print_vector_round(ids, stdout);
    err = igraph_is_potentially_connected_degree_sequence(ods, ids, &res);
    if (err != IGRAPH_SUCCESS) {
        printf("error!\n"); goto cleanup;
    }
    printf("%s\n", res ? "true" : "false");

cleanup:
    igraph_vector_destroy(ods);
    igraph_vector_destroy(ids);
}

int main() {
    igraph_vector_t ds, ods, ids;

    printf("Undirected:\n");

    /* Null graph: not potentially connected. */
    igraph_vector_init_int_end(&ds, -1, -1);
    pc_undirected_print_destroy(&ds);

    igraph_vector_init_int_end(&ds, -1, 0, -1);
    pc_undirected_print_destroy(&ds);

    igraph_vector_init_int_end(&ds, -1, 1, -1);
    pc_undirected_print_destroy(&ds);

    /* Not potentially connected. */
    igraph_vector_init_int_end(&ds, -1, 1, 1, 1, 1, -1);
    pc_undirected_print_destroy(&ds);

    igraph_vector_init_int_end(&ds, -1, 1, 2, 3, -1);
    pc_undirected_print_destroy(&ds);

    igraph_vector_init_int_end(&ds, -1, 1, 2, 3, 0, -1);
    pc_undirected_print_destroy(&ds);

    /* Non-even sum */
    igraph_vector_init_int_end(&ds, -1, 3, 2, -1);
    pc_undirected_print_destroy(&ds);

    /* Negative value */
    igraph_vector_init_int_end(&ds, -1, -2, 2, -1);
    pc_undirected_print_destroy(&ds);

    printf("\n\nDirected:\n");

    igraph_vector_init_int_end(&ods, -1, -1);
    igraph_vector_init_int_end(&ids, -1, -1);
    pc_directed_print_destroy(&ods, &ids);

    igraph_vector_init_int_end(&ods, -1, 0, -1);
    igraph_vector_init_int_end(&ids, -1, 0, -1);
    pc_directed_print_destroy(&ods, &ids);

    igraph_vector_init_int_end(&ods, -1, 0, -1);
    igraph_vector_init_int_end(&ids, -1, 1, -1);
    pc_directed_print_destroy(&ods, &ids);

    igraph_vector_init_int_end(&ods, -1, 1, 2, 3, -1);
    igraph_vector_init_int_end(&ids, -1, 2, 2, 2, -1);
    pc_directed_print_destroy(&ods, &ids);

    igraph_vector_init_int_end(&ods, -1, 1, 2, 3, -1);
    igraph_vector_init_int_end(&ids, -1, 2, 4, 0, -1);
    pc_directed_print_destroy(&ods, &ids);

    igraph_vector_init_int_end(&ods, -1, 1, 2, 3, -1);
    igraph_vector_init_int_end(&ids, -1, 2, 6, -2, -1);
    pc_directed_print_destroy(&ods, &ids);

    VERIFY_FINALLY_STACK();

    return 0;
}
