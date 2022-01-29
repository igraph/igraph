#include <igraph.h>

int main() {

    igraph_t g;

    igraph_attribute_combination_t comb;

    igraph_set_attribute_table(&igraph_cattribute_table);

    igraph_small(&g, 2, IGRAPH_DIRECTED,
                 0, 1, 0, 1, -1);

    SETEAB(&g, "type", 0, 1);
    SETEAB(&g, "type", 1, 0);

    igraph_attribute_combination(&comb,
                                 "weight", IGRAPH_ATTRIBUTE_COMBINE_SUM,
                                 "type",   IGRAPH_ATTRIBUTE_COMBINE_FIRST,
                                 "",       IGRAPH_ATTRIBUTE_COMBINE_IGNORE,
                                 IGRAPH_NO_MORE_ATTRIBUTES);
    igraph_simplify(&g, /*multiple=*/ 1, /*loops=*/ 1, &comb);
    igraph_write_graph_graphml(&g, stdout, /*prefixattr=*/ 1);
    igraph_destroy(&g);

    return 0;
}
