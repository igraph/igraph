/* vim:set ts=4 sw=4 sts=4 et: */

#include <igraph.h>

#include "test_utilities.h"

int main(void) {
    igraph_t g;
    igraph_vector_t result;
    igraph_vector_t weights;

    igraph_vector_init(&result, 0);

    /* null graph */
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    igraph_vector_init(&weights, 0);

    printf("Null graph:\n");
    igraph_diversity(&g, &weights, &result, igraph_vss_all());
    print_vector(&result);

    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    /* graph with no edges */
    igraph_empty(&g, 5, IGRAPH_UNDIRECTED);
    igraph_vector_init(&weights, 0);

    printf("Empty graph:\n");
    igraph_diversity(&g, &weights, &result, igraph_vss_all());
    print_vector(&result);

    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    /* real graph */
    igraph_small(&g, 4, IGRAPH_UNDIRECTED, 0,1, 0,2, 1,2, 1,3, 2,3, -1);
    igraph_vector_init_int_end(&weights, -1, 3, 2, 8, 1, 1, -1);

    printf("Graph with 4 nodes and 5 edges:\n");
    igraph_diversity(&g, &weights, &result, igraph_vss_all());
    print_vector(&result);

    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    /* degree-one vertices */
    igraph_kary_tree(&g, 10, 2, IGRAPH_TREE_UNDIRECTED);
    igraph_vector_init_range(&weights, 1, igraph_ecount(&g) + 1);

    printf("Tree (having degree-one vertices):\n");
    igraph_diversity(&g, &weights, &result, igraph_vss_all());
    print_vector(&result);

    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    /* error conditions are tested from now on */
    VERIFY_FINALLY_STACK();

    /* graph with multiple edges */
    igraph_small(&g, 3, IGRAPH_UNDIRECTED, 0,1, 0,2, 2,0, -1);
    igraph_vector_init_int_end(&weights, -1, 3, 2, 8, -1);
    CHECK_ERROR(igraph_diversity(&g, &weights, &result, igraph_vss_all()), IGRAPH_EINVAL);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    /* negative weights */
    igraph_small(&g, 3, IGRAPH_UNDIRECTED, 0,1, 0,2, 2,1, -1);
    igraph_vector_init_int_end(&weights, -1, 3, -2, 8, -1);
    CHECK_ERROR(igraph_diversity(&g, &weights, &result, igraph_vss_all()), IGRAPH_EINVAL);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    /* directed graph */
    igraph_small(&g, 3, IGRAPH_DIRECTED, 0,1, 0,2, -1);
    igraph_vector_init_int_end(&weights, -1, 3, 2, -1);
    CHECK_ERROR(igraph_diversity(&g, &weights, &result, igraph_vss_all()), IGRAPH_EINVAL);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    igraph_vector_destroy(&result);

    VERIFY_FINALLY_STACK();

    return 0;
}
