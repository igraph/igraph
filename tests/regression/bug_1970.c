
#include <igraph.h>

int main(void) {
    FILE *file;
    igraph_t graph;
    igraph_error_t err;

    igraph_set_attribute_table(&igraph_cattribute_table);

    igraph_set_error_handler(igraph_error_handler_printignore);

    file = fopen("bug_1970.graphml", "r");
    IGRAPH_ASSERT(file != NULL);

    err = igraph_read_graph_graphml(&graph, file, 0);
    fclose(file);
    IGRAPH_ASSERT(err != IGRAPH_SUCCESS);

    /* graph creation should have failed, no need to destroy 'graph' */

    return 0;
}
