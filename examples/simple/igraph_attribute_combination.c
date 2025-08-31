#include <igraph.h>

int main(void) {
    igraph_t graph;
    igraph_attribute_combination_t comb;

    /* Initialize the library. */
    igraph_setup();

    igraph_set_attribute_table(&igraph_cattribute_table);

    igraph_small(&graph, 2, IGRAPH_DIRECTED,
                 0,1, 0,1,
                 -1);

    SETEAB(&graph, "type", 0, true);
    SETEAB(&graph, "type", 1, false);

    igraph_attribute_combination(&comb,
                                 "weight", IGRAPH_ATTRIBUTE_COMBINE_SUM,
                                 "type",   IGRAPH_ATTRIBUTE_COMBINE_FIRST,
                                 "",       IGRAPH_ATTRIBUTE_COMBINE_IGNORE,
                                 IGRAPH_NO_MORE_ATTRIBUTES);
    igraph_simplify(&graph, /*remove_multiple=*/ true, /*remove_loops=*/ true, &comb);
    igraph_write_graph_graphml(&graph, stdout, /*prefixattr=*/ true);

    igraph_destroy(&graph);
    igraph_attribute_combination_destroy(&comb);

    return 0;
}
