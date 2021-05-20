#include <igraph.h>

#include "../unit/test_utilities.inc"

/* Regression test for https://github.com/igraph/igraph/issues/1760 */

int test_unweighted(const igraph_t* g, igraph_integer_t from, const igraph_vs_t* to) {
    igraph_vector_ptr_t vpath, epath;
    igraph_integer_t i;
    igraph_integer_t num_paths;
    igraph_vector_long_t predecessors;
    igraph_vector_long_t inbound_edges;

    printf("Unweighted case\n");
    printf("---------------\n\n");

    IGRAPH_CHECK(igraph_vs_size(g, to, &num_paths));
    IGRAPH_CHECK(igraph_vector_ptr_init(&vpath, num_paths));
    IGRAPH_CHECK(igraph_vector_ptr_init(&epath, num_paths));
    IGRAPH_CHECK(igraph_vector_long_init(&predecessors, 0));
    IGRAPH_CHECK(igraph_vector_long_init(&inbound_edges, 0));

    for (i = 0; i < igraph_vector_ptr_size(&vpath); i++) {
        VECTOR(vpath)[i] = igraph_Calloc(1, igraph_vector_t);
        VECTOR(epath)[i] = igraph_Calloc(1, igraph_vector_t);
        IGRAPH_CHECK(igraph_vector_init(VECTOR(vpath)[i], 0));
        IGRAPH_CHECK(igraph_vector_init(VECTOR(epath)[i], 0));
    }

    IGRAPH_CHECK(igraph_get_shortest_paths(
        g, &vpath, &epath, from, *to, IGRAPH_IN,
        &predecessors, &inbound_edges
    ));

    printf("Vertices:\n");
    for (i = 0; i < igraph_vector_ptr_size(&vpath); i++) {
        print_vector(VECTOR(vpath)[i]);
        igraph_vector_destroy(VECTOR(vpath)[i]);
    }
    printf("\n");

    printf("Edges:\n");
    for (i = 0; i < igraph_vector_ptr_size(&epath); i++) {
        print_vector(VECTOR(epath)[i]);
        igraph_vector_destroy(VECTOR(epath)[i]);
    }
    printf("\n");

    printf("Predecessors:\n");
    print_vector_long(&predecessors);
    printf("\n");

    printf("Inbound edges:\n");
    print_vector_long(&inbound_edges);
    printf("\n");

    igraph_vector_long_destroy(&inbound_edges);
    igraph_vector_long_destroy(&predecessors);
    igraph_vector_ptr_destroy_all(&epath);
    igraph_vector_ptr_destroy_all(&vpath);

    return IGRAPH_SUCCESS;
}

int test_weighted(
    const igraph_t* g, const igraph_vector_t* weights, igraph_integer_t from,
    const igraph_vs_t* to, igraph_bool_t use_bellman_ford
) {
    igraph_vector_ptr_t vpath, epath;
    igraph_integer_t i;
    igraph_integer_t num_paths;
    igraph_vector_long_t predecessors;
    igraph_vector_long_t inbound_edges;

    printf("Weighted case\n");
    printf("-------------\n\n");

    printf("Algorithm: %s\n\n", use_bellman_ford ? "Bellman-Ford" : "Dijkstra");

    IGRAPH_CHECK(igraph_vs_size(g, to, &num_paths));
    IGRAPH_CHECK(igraph_vector_ptr_init(&vpath, num_paths));
    IGRAPH_CHECK(igraph_vector_ptr_init(&epath, num_paths));
    IGRAPH_CHECK(igraph_vector_long_init(&predecessors, 0));
    IGRAPH_CHECK(igraph_vector_long_init(&inbound_edges, 0));

    for (i = 0; i < igraph_vector_ptr_size(&vpath); i++) {
        VECTOR(vpath)[i] = igraph_Calloc(1, igraph_vector_t);
        VECTOR(epath)[i] = igraph_Calloc(1, igraph_vector_t);
        IGRAPH_CHECK(igraph_vector_init(VECTOR(vpath)[i], 0));
        IGRAPH_CHECK(igraph_vector_init(VECTOR(epath)[i], 0));
    }

    if (use_bellman_ford) {
        IGRAPH_CHECK(igraph_get_shortest_paths_bellman_ford(
            g, &vpath, &epath, from, *to, weights, IGRAPH_IN,
            &predecessors, &inbound_edges
        ));
    } else {
        IGRAPH_CHECK(igraph_get_shortest_paths_dijkstra(
            g, &vpath, &epath, from, *to, weights, IGRAPH_IN,
            &predecessors, &inbound_edges
        ));
    }

    printf("Vertices:\n");
    for (i = 0; i < igraph_vector_ptr_size(&vpath); i++) {
        print_vector(VECTOR(vpath)[i]);
        igraph_vector_destroy(VECTOR(vpath)[i]);
    }
    printf("\n");

    printf("Edges:\n");
    for (i = 0; i < igraph_vector_ptr_size(&epath); i++) {
        print_vector(VECTOR(epath)[i]);
        igraph_vector_destroy(VECTOR(epath)[i]);
    }
    printf("\n");

    printf("Predecessors:\n");
    print_vector_long(&predecessors);
    printf("\n");

    printf("Inbound edges:\n");
    print_vector_long(&inbound_edges);
    printf("\n");

    igraph_vector_long_destroy(&inbound_edges);
    igraph_vector_long_destroy(&predecessors);
    igraph_vector_ptr_destroy_all(&epath);
    igraph_vector_ptr_destroy_all(&vpath);

    return IGRAPH_SUCCESS;
}

int main(int argc, char* argv[]) {
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
