
#include <igraph.h>
#include <stdio.h>

#include "test_utilities.h"

void test_input(const char *filename) {
    igraph_t graph;
    FILE *ifile;

    ifile = fopen(filename, "r");
    if (! ifile) {
        fprintf(stderr, "Cannot open '%s'.\n", filename);
        abort();
    }

    printf("===== %s =====\n", filename);

    IGRAPH_ASSERT(igraph_read_graph_gml(&graph, ifile) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_write_graph_gml(&graph, stdout, NULL, "") == IGRAPH_SUCCESS);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    printf("======================\n\n");
}

int main() {

    /* Enable attribute handling */
    igraph_set_attribute_table(&igraph_cattribute_table);

    test_input("graph1.gml");
    test_input("graph2.gml");
    test_input("graph3.gml");

    return 0;
}
