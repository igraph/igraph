
#include <igraph.h>

#include "test_utilities.inc"

#define BIGRAPHICAL_PRINT_DESTROY(deg1, deg2) \
    igraph_is_bigraphical(&(deg1), &(deg2), IGRAPH_SIMPLE_SW, &simple); \
    igraph_is_bigraphical(&(deg1), &(deg2), IGRAPH_MULTI_SW, &multi); \
    print_vector_round(&(deg1)); \
    print_vector_round(&(deg2)); \
    printf("simple: %s, multi: %s\n\n", simple ? "true" : "false", multi ? "true" : "false"); \
    igraph_vector_destroy(&(deg1)); \
    igraph_vector_destroy(&(deg2));

int main() {
    igraph_vector_t deg1, deg2;
    igraph_bool_t simple, multi;

    igraph_vector_init(&deg1, 0);
    igraph_vector_init(&deg2, 0);
    BIGRAPHICAL_PRINT_DESTROY(deg1, deg2);

    igraph_vector_init_int_end(&deg1, -1, 3, 3, -1);
    igraph_vector_init_int_end(&deg2, -1, 1, 2, 3, -1);
    BIGRAPHICAL_PRINT_DESTROY(deg1, deg2);

    igraph_vector_init_int_end(&deg1, -1, 3, 2, 1, -1);
    igraph_vector_init_int_end(&deg2, -1, 1, 2, 3, -1);
    BIGRAPHICAL_PRINT_DESTROY(deg1, deg2);

    igraph_vector_init_int_end(&deg1, -1, 1, 1, 1, 1, -1);
    igraph_vector_init_int_end(&deg2, -1, 2, 3, -1);
    BIGRAPHICAL_PRINT_DESTROY(deg1, deg2);

    igraph_vector_init_int_end(&deg1, -1, 1, 1, 1, 1, -1);
    igraph_vector_init_int_end(&deg2, -1, 2, 2, -1);
    BIGRAPHICAL_PRINT_DESTROY(deg1, deg2);

    igraph_vector_init_int_end(&deg1, -1, 1, 2, 0, 3, 0, -1);
    igraph_vector_init_int_end(&deg2, -1, 2, 3, 1, -1);
    BIGRAPHICAL_PRINT_DESTROY(deg1, deg2);

    VERIFY_FINALLY_STACK();

    return 0;
}
