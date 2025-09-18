
#include <igraph.h>
#include <stdlib.h>

int main(void) {

    igraph_t ring, g, *component;
    igraph_graph_list_t complist;
    igraph_int_t i;

    /* Initialize the library. */
    igraph_setup();

    igraph_graph_list_init(&complist, 0);

    /* A ring, a single component */
    igraph_ring(&ring, 10, IGRAPH_UNDIRECTED, 0, 1);

    igraph_decompose(&ring, &complist, IGRAPH_WEAK, -1, 0);
    component = igraph_graph_list_get_ptr(&complist, 0);
    igraph_write_graph_edgelist(component, stdout);
    igraph_destroy(&ring);
    igraph_graph_list_clear(&complist);

    /* Random graph with a giant component */
    igraph_erdos_renyi_game_gnp(&g, 100, 4.0 / 100, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    igraph_decompose(&g, &complist, IGRAPH_WEAK, -1, 20);
    if (igraph_graph_list_size(&complist) != 1) {
        return 1;
    }
    igraph_destroy(&g);
    igraph_graph_list_clear(&complist);

    /* A toy graph, three components maximum, with at least 2 vertices each */
    igraph_small(&g, 0, IGRAPH_DIRECTED,
                 0, 1, 1, 2, 2, 0,
                 3, 4, 4, 5, 5, 6,
                 8, 9, 9, 10,
                 -1);
    igraph_decompose(&g, &complist, IGRAPH_WEAK, 3, 2);
    for (i = 0; i < igraph_graph_list_size(&complist); i++) {
        component = igraph_graph_list_get_ptr(&complist, i);
        igraph_write_graph_edgelist(component, stdout);
    }
    igraph_destroy(&g);

    igraph_graph_list_destroy(&complist);

    return 0;
}
