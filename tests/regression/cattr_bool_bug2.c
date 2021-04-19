
#include <igraph.h>
#include <stdio.h>

#include "../unit/test_utilities.inc"

#define FILENAME "mybool.graphml.xml"

int main() {

    igraph_t graph;
    igraph_error_handler_t* oldhandler;
    int result;
    FILE* ifile = fopen("cattr_bool_bug2.graphml", "r");

    if (!ifile) {
        printf("Cannot open input file\n");
        return 1;
    }

    igraph_set_attribute_table(&igraph_cattribute_table);

    oldhandler = igraph_set_error_handler(igraph_error_handler_ignore);
    if ((result = igraph_read_graph_graphml(&graph, ifile, 0))) {
        /* maybe it is simply disabled at compile-time */
        if (result == IGRAPH_UNIMPLEMENTED) {
            return 77;
        }
        printf("Failed to read GraphML file\n");
        return 1;
    }
    igraph_set_error_handler(oldhandler);

    fclose(ifile);

    IGRAPH_ASSERT(igraph_cattribute_has_attr(&graph, IGRAPH_ATTRIBUTE_GRAPH, "mybool"));

    /* Boolean attribute value is expected to be true */
    IGRAPH_ASSERT(igraph_cattribute_GAB(&graph, "mybool"));

    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
