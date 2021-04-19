#include <igraph.h>

int main() {

    igraph_t g;
    igraph_strvector_t names, weights;
    int i;
    char str[2] = " ";
    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_small(&g, 7, IGRAPH_UNDIRECTED, 0, 1, 0, 2, 1, 2, 1, 3, 2, 4, 3, 4, -1);


    printf("Output without isolates:\n");
    igraph_write_graph_lgl(&g, stdout, /*names*/ NULL, /*weights*/ NULL, /*isolates*/ 0);


    printf("\nOutput with isolates:\n");
    igraph_write_graph_lgl(&g, stdout, /*names*/ NULL, /*weights*/ NULL, /*isolates*/ 1);


    printf("\nOutput vertex and edge labels:\n");
    igraph_strvector_init(&names, 7);
    for (i = 0; i < 7; i++) {
        str[0] = 'A' + i;
        igraph_strvector_set(&names, i, str);
    }
    SETVASV(&g, "names", &names);

    igraph_strvector_init(&weights, 6);
    for (i = 0; i < 6; i++) {
        str[0] = '3' + i;
        igraph_strvector_set(&weights, i, str);
    }
    SETEASV(&g, "weights", &weights);

    igraph_write_graph_lgl(&g, stdout, "names", "weights", /*isolates*/ 0);

    igraph_strvector_destroy(&names);
    igraph_strvector_destroy(&weights);
    igraph_destroy(&g);

    return 0;
}
