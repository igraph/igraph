#include <igraph.h>

int main(void) {

    igraph_t graph;
    igraph_strvector_t names;
    igraph_vector_t weights;
    igraph_integer_t i;
    igraph_integer_t vcount, ecount;

    igraph_set_attribute_table(&igraph_cattribute_table);

    igraph_small(&graph, 7, IGRAPH_UNDIRECTED,
                 0,1, 1,3, 1,2, 2,0, 4,2, 3,4,
                 -1);
    vcount = igraph_vcount(&graph);
    ecount = igraph_ecount(&graph);


    printf("Output without isolates:\n");
    igraph_write_graph_lgl(&graph, stdout, /*names*/ NULL, /*weights*/ NULL, /*isolates*/ 0);


    printf("\nOutput with isolates:\n");
    igraph_write_graph_lgl(&graph, stdout, /*names*/ NULL, /*weights*/ NULL, /*isolates*/ 1);


    printf("\nOutput vertex and edge labels:\n");
    igraph_strvector_init(&names, vcount);
    for (i = 0; i < vcount; i++) {
        char str[2] = " "; /* initialize to ensure presence of null terminator */
        str[0] = 'A' + i;
        igraph_strvector_set(&names, i, str);
    }
    SETVASV(&graph, "names", &names);

    igraph_vector_init_range(&weights, 1, ecount + 1);
    SETEANV(&graph, "weights", &weights);

    igraph_write_graph_lgl(&graph, stdout, "names", "weights", /*isolates*/ 0);

    igraph_strvector_destroy(&names);
    igraph_vector_destroy(&weights);
    igraph_destroy(&graph);

    return 0;
}
