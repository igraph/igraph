
#include <igraph.h>
#include <stdio.h>

#define FILENAME "mybool.graphml.xml"

int main() {

    igraph_t graph;
    igraph_error_handler_t* oldhandler;
    int result;
    FILE* ifile = fopen("cattr_bool_bug2.graphml", "r");

    if (!ifile) {
        printf("Cannot open input file");
        return 1;
    }

    igraph_i_set_attribute_table(&igraph_cattribute_table);

    oldhandler = igraph_set_error_handler(igraph_error_handler_ignore);
    if ((result = igraph_read_graph_graphml(&graph, ifile, 0))) {
        /* maybe it is simply disabled at compile-time */
        if (result == IGRAPH_UNIMPLEMENTED) {
            return 77;
        }
        return 1;
    }
    igraph_set_error_handler(oldhandler);

    fclose(ifile);

    if (!igraph_cattribute_has_attr(&graph, IGRAPH_ATTRIBUTE_GRAPH, "mybool")) {
        printf("boolean value mybool not found\n");
        return 2;
    } else {
        igraph_bool_t value = igraph_cattribute_GAB(&graph, "mybool");
        printf("found boolean value %d\n", value);
    }

    igraph_destroy(&graph);

    return 0;
}
