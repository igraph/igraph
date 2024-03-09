
#include <igraph.h>

int main(void) {
    FILE *file;
    igraph_t graph;

    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_set_error_handler(igraph_error_handler_ignore);
    igraph_set_warning_handler(igraph_warning_handler_ignore);

    file = fopen("bug_2497.gml", "r");
    IGRAPH_ASSERT(file != NULL);

    IGRAPH_ASSERT(IGRAPH_FINALLY_STACK_EMPTY);
    IGRAPH_ASSERT(igraph_read_graph_gml(&graph, file) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(IGRAPH_FINALLY_STACK_EMPTY);

    fclose(file);

    igraph_write_graph_gml(&graph, stdout, IGRAPH_WRITE_GML_DEFAULT_SW, NULL, "");
    IGRAPH_ASSERT(IGRAPH_FINALLY_STACK_EMPTY);

    igraph_destroy(&graph);

    return 0;
}
