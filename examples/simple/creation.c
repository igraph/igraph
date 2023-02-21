
#include <igraph.h>
#include <assert.h>

int main(void) {
    igraph_t graph;
    igraph_vector_int_t edges;

    /* Create a directed graph with no vertices or edges. */
    igraph_empty(&graph, 0, IGRAPH_DIRECTED);

    /* Add 5 vertices. Vertex IDs will range from 0 to 4, inclusive. */
    igraph_add_vertices(&graph, 5, NULL);

    /* Add 5 edges, specified as 5 consecutive pairs of vertex IDs
     * stored in an integer vector. */
    igraph_vector_int_init_int(&edges, 10,
                               0,1, 0,2, 3,1, 2,1, 0,4);
    igraph_add_edges(&graph, &edges, NULL);

    igraph_vector_int_destroy(&edges);

    /* Now the graph has 5 vertices and 5 edges. */
    assert(igraph_vcount(&graph) == 5);
    assert(igraph_ecount(&graph) == 5);

    igraph_destroy(&graph);

    return 0;
}
