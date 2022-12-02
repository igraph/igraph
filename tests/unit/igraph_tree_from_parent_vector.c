
#include <igraph.h>
#include "test_utilities.h"

int main(void) {
    igraph_t graph;
    igraph_vector_int_t parents;
    igraph_bool_t is_tree, is_forest;

    igraph_vector_int_init_int(&parents, 5,
                               4, 4, 1, -2, 3);

    printf("Out-tree:\n");
    igraph_tree_from_parent_vector(&graph, &parents, IGRAPH_TREE_OUT);
    print_graph(&graph);
    igraph_is_tree(&graph, &is_tree, NULL, IGRAPH_OUT);
    IGRAPH_ASSERT(is_tree);
    igraph_destroy(&graph);

    printf("\nIn-tree:\n");
    igraph_tree_from_parent_vector(&graph, &parents, IGRAPH_TREE_IN);
    print_graph(&graph);
    igraph_is_tree(&graph, &is_tree, NULL, IGRAPH_IN);
    IGRAPH_ASSERT(is_tree);
    igraph_destroy(&graph);

    printf("\nUndirected tree:\n");
    igraph_tree_from_parent_vector(&graph, &parents, IGRAPH_TREE_UNDIRECTED);
    print_graph(&graph);
    igraph_is_tree(&graph, &is_tree, NULL, IGRAPH_ALL);
    IGRAPH_ASSERT(is_tree);
    igraph_destroy(&graph);

    printf("\nForest:\n");
    VECTOR(parents)[0] = -1;
    igraph_tree_from_parent_vector(&graph, &parents, IGRAPH_TREE_OUT);
    print_graph(&graph);
    igraph_is_forest(&graph, &is_forest, NULL, IGRAPH_OUT);
    IGRAPH_ASSERT(is_forest);
    igraph_destroy(&graph);

    /* Invalid tree mode */
    CHECK_ERROR(igraph_tree_from_parent_vector(&graph, &parents, -1), IGRAPH_EINVAL);

    /* Invalid vertex */
    VECTOR(parents)[0] = 5;
    CHECK_ERROR(igraph_tree_from_parent_vector(&graph, &parents, IGRAPH_TREE_OUT), IGRAPH_EINVVID);

    /* Self-loop */
    VECTOR(parents)[0] = 0;
    CHECK_ERROR(igraph_tree_from_parent_vector(&graph, &parents, IGRAPH_TREE_OUT), IGRAPH_EINVAL);

    /* Longer cycle */
    VECTOR(parents)[0] = 4;
    VECTOR(parents)[3] = 0;
    CHECK_ERROR(igraph_tree_from_parent_vector(&graph, &parents, IGRAPH_TREE_OUT), IGRAPH_EINVAL);

    printf("\nEdgeless graph:\n");
    igraph_vector_int_fill(&parents, -1);
    igraph_tree_from_parent_vector(&graph, &parents, IGRAPH_TREE_OUT);
    print_graph(&graph);
    IGRAPH_ASSERT(igraph_ecount(&graph) == 0);
    igraph_destroy(&graph);

    printf("\nNull graph:\n");
    igraph_vector_int_clear(&parents);
    igraph_tree_from_parent_vector(&graph, &parents, IGRAPH_TREE_OUT);
    print_graph(&graph);
    IGRAPH_ASSERT(igraph_vcount(&graph) == 0);
    igraph_destroy(&graph);

    igraph_vector_int_destroy(&parents);

    VERIFY_FINALLY_STACK();

    return 0;
}
