
#include <igraph.h>
#include <stdio.h>
#include <stdlib.h>

void check_attr(igraph_t *graph, int offset) {

    if (!igraph_cattribute_has_attr(graph, IGRAPH_ATTRIBUTE_GRAPH, "name")) {
        printf("No graph attribute `name`\n");
        exit(offset + 2);
    }
    if (!igraph_cattribute_has_attr(graph, IGRAPH_ATTRIBUTE_GRAPH, "type")) {
        printf("No graph attribute `type`\n");
        exit(offset + 3);
    }
    if (!igraph_cattribute_has_attr(graph, IGRAPH_ATTRIBUTE_GRAPH, "p")) {
        printf("No graph attribute `p`\n");
        exit(offset + 4);
    }
    if (!igraph_cattribute_has_attr(graph, IGRAPH_ATTRIBUTE_VERTEX,
                                    "name")) {
        printf("No vertex attribute `id`\n");
        exit(offset + 5);
    }
    if (!igraph_cattribute_has_attr(graph, IGRAPH_ATTRIBUTE_EDGE,
                                    "weight")) {
        printf("No edge attribute `weight'\n");
        exit(offset + 6);
    }
}

int main() {

    igraph_t graph;
    igraph_error_handler_t* oldhandler;
    int result;
    FILE *ifile = fopen("cattr_bool_bug.graphml", "r");

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

    check_attr(&graph, 10);

    igraph_to_directed(&graph, IGRAPH_TO_DIRECTED_ARBITRARY);

    check_attr(&graph, 20);

    if (GAB(&graph, "loops")) {
        return 2;
    }

    igraph_destroy(&graph);

    return 0;
}
