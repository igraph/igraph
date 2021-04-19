
#include <igraph.h>
#include <stdlib.h>

void free_complist(igraph_vector_ptr_t *complist) {
    long int i;
    for (i = 0; i < igraph_vector_ptr_size(complist); i++) {
        igraph_destroy(VECTOR(*complist)[i]);
        igraph_free(VECTOR(*complist)[i]);
    }
}

int main() {

    igraph_t ring, g;
    igraph_vector_ptr_t complist;
    long int i;
    igraph_real_t edges[] = { 0, 1, 1, 2, 2, 0,
                              3, 4, 4, 5, 5, 6,
                              8, 9, 9, 10
                            };
    igraph_vector_t v;

    /* A ring, a single component */
    igraph_ring(&ring, 10, IGRAPH_UNDIRECTED, 0, 1);

    igraph_vector_ptr_init(&complist, 0);
    igraph_decompose(&ring, &complist, IGRAPH_WEAK, -1, 0);
    igraph_write_graph_edgelist(VECTOR(complist)[0], stdout);
    free_complist(&complist);
    igraph_destroy(&ring);

    /* Random graph with a giant component */
    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 100, 4.0 / 100,
                            IGRAPH_UNDIRECTED, 0);
    igraph_decompose(&g, &complist, IGRAPH_WEAK, -1, 20);
    if (igraph_vector_ptr_size(&complist) != 1) {
        return 1;
    }
    free_complist(&complist);
    igraph_destroy(&g);

    /* A toy graph, three components maximum, with at least 2 vertices each */
    igraph_create(&g,
                  igraph_vector_view(&v, edges, sizeof(edges) / sizeof(igraph_real_t)),
                  0, IGRAPH_DIRECTED);
    igraph_decompose(&g, &complist, IGRAPH_WEAK, 3, 2);
    for (i = 0; i < igraph_vector_ptr_size(&complist); i++) {
        igraph_write_graph_edgelist(VECTOR(complist)[i], stdout);
    }
    free_complist(&complist);
    igraph_destroy(&g);

    igraph_vector_ptr_destroy(&complist);

    return 0;
}
