#include <igraph.h>

#include "../unit/test_utilities.h"

/* Regression test for https://github.com/igraph/igraph/issues/1814 */

void test_igraph_to_undirected(igraph_to_undirected_t mode) {
    igraph_t g;
    igraph_attribute_combination_t comb;

    igraph_small(&g, 4, IGRAPH_DIRECTED,
                 0,1, 0,1, 1,0,
                 1,1, 1,1,
                 -1);

    SETEAN(&g, "weight", 0, 1);
    SETEAN(&g, "weight", 1, 2);
    SETEAN(&g, "weight", 2, 3);
    SETEAN(&g, "weight", 3, 4);
    SETEAN(&g, "weight", 4, 2);

    igraph_attribute_combination(
        &comb, "weight", IGRAPH_ATTRIBUTE_COMBINE_SUM,
        IGRAPH_NO_MORE_ATTRIBUTES
    );
    igraph_to_undirected(&g, mode, &comb);
    igraph_attribute_combination_destroy(&comb);

    igraph_write_graph_gml(&g, stdout, IGRAPH_WRITE_GML_DEFAULT_SW, 0, "");

    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
}

int main(void) {
    igraph_set_attribute_table(&igraph_cattribute_table);

    printf("to_undirected(COLLAPSE)\n");
    printf("=======================\n\n");
    test_igraph_to_undirected(IGRAPH_TO_UNDIRECTED_COLLAPSE);
    printf("\n");

    printf("to_undirected(MUTUAL)\n");
    printf("=====================\n\n");
    test_igraph_to_undirected(IGRAPH_TO_UNDIRECTED_MUTUAL);
    printf("\n");

    printf("to_undirected(EACH)\n");
    printf("===================\n\n");
    test_igraph_to_undirected(IGRAPH_TO_UNDIRECTED_EACH);
    printf("\n");

    return 0;
}
