#include <igraph.h>

#include "../unit/test_utilities.h"

/* Regression test for https://github.com/igraph/igraph/issues/1760 */

int test_unweighted(const igraph_t* g, igraph_integer_t from, const igraph_vs_t* to) {
    igraph_vector_int_list_t vpath, epath;
    igraph_integer_t num_paths;
    igraph_vector_int_t parents;
    igraph_vector_int_t inbound_edges;

    printf("Unweighted case\n");
    printf("---------------\n\n");

    IGRAPH_CHECK(igraph_vs_size(g, to, &num_paths));
    IGRAPH_CHECK(igraph_vector_int_list_init(&vpath, 0));
    IGRAPH_CHECK(igraph_vector_int_list_init(&epath, 0));
    IGRAPH_CHECK(igraph_vector_int_init(&parents, 0));
    IGRAPH_CHECK(igraph_vector_int_init(&inbound_edges, 0));

    IGRAPH_CHECK(igraph_get_shortest_paths(
        g, &vpath, &epath, from, *to, IGRAPH_IN,
        &parents, &inbound_edges
    ));

    printf("Vertices:\n");
    print_vector_int_list(&vpath);
    printf("\n");

    printf("Edges:\n");
    print_vector_int_list(&epath);
    printf("\n");

    printf("Parents:\n");
    print_vector_int(&parents);
    printf("\n");

    printf("Inbound edges:\n");
    print_vector_int(&inbound_edges);
    printf("\n");

    igraph_vector_int_destroy(&inbound_edges);
    igraph_vector_int_destroy(&parents);
    igraph_vector_int_list_destroy(&epath);
    igraph_vector_int_list_destroy(&vpath);

    return IGRAPH_SUCCESS;
}

int test_weighted(
    const igraph_t* g, const igraph_vector_t* weights, igraph_integer_t from,
    const igraph_vs_t* to, igraph_bool_t use_bellman_ford
) {
    igraph_vector_int_list_t vpath, epath;
    igraph_integer_t num_paths;
    igraph_vector_int_t parents;
    igraph_vector_int_t inbound_edges;

    printf("Weighted case\n");
    printf("-------------\n\n");

    printf("Algorithm: %s\n\n", use_bellman_ford ? "Bellman-Ford" : "Dijkstra");

    IGRAPH_CHECK(igraph_vs_size(g, to, &num_paths));
    IGRAPH_CHECK(igraph_vector_int_list_init(&vpath, 0));
    IGRAPH_CHECK(igraph_vector_int_list_init(&epath, 0));
    IGRAPH_CHECK(igraph_vector_int_init(&parents, 0));
    IGRAPH_CHECK(igraph_vector_int_init(&inbound_edges, 0));

    if (use_bellman_ford) {
        IGRAPH_CHECK(igraph_get_shortest_paths_bellman_ford(
            g, &vpath, &epath, from, *to, weights, IGRAPH_IN,
            &parents, &inbound_edges
        ));
    } else {
        IGRAPH_CHECK(igraph_get_shortest_paths_dijkstra(
            g, &vpath, &epath, from, *to, weights, IGRAPH_IN,
            &parents, &inbound_edges
        ));
    }

    printf("Vertices:\n");
    print_vector_int_list(&vpath);
    printf("\n");

    printf("Edges:\n");
    print_vector_int_list(&epath);
    printf("\n");

    printf("Parents:\n");
    print_vector_int(&parents);
    printf("\n");

    printf("Inbound edges:\n");
    print_vector_int(&inbound_edges);
    printf("\n");

    igraph_vector_int_destroy(&inbound_edges);
    igraph_vector_int_destroy(&parents);
    igraph_vector_int_list_destroy(&epath);
    igraph_vector_int_list_destroy(&vpath);

    return IGRAPH_SUCCESS;
}

int main(void) {
    igraph_t g;
    igraph_vector_t weights;
    igraph_vs_t to;

    igraph_set_warning_handler(igraph_warning_handler_ignore);

    igraph_small(&g, 4, /* directed = */ 1, 0, 1, 1, 2, 1, 3, -1);
    igraph_vs_vector_small(&to, 0, 3, -1);
    igraph_vector_init(&weights, 3);
    VECTOR(weights)[0] = 1;
    VECTOR(weights)[1] = 2;
    VECTOR(weights)[2] = 2;

    /* Test unweighted case */
    if (test_unweighted(&g, 2, &to)) {
        return 1;
    }

    /* Test weighted case */
    if (test_weighted(&g, &weights, 2, &to, /* use_bellman_ford = */ 0)) {
        return 2;
    }
    if (test_weighted(&g, &weights, 2, &to, /* use_bellman_ford = */ 1)) {
        return 3;
    }

    igraph_vector_destroy(&weights);
    igraph_vs_destroy(&to);
    igraph_destroy(&g);

    return 0;
}
