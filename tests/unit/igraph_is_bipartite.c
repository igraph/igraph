
#include <igraph.h>

#include "test_utilities.inc"

int main() {

    igraph_t graph;
    igraph_bool_t bipartite;
    igraph_vector_bool_t types;

    igraph_vector_bool_init(&types, 5);

    /* Null graph */
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    igraph_is_bipartite(&graph, &bipartite, &types);
    IGRAPH_ASSERT(bipartite);
    IGRAPH_ASSERT(igraph_vector_bool_size(&types) == igraph_vcount(&graph));
    igraph_destroy(&graph);

    /* Singleton graph */
    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    igraph_is_bipartite(&graph, &bipartite, &types);
    IGRAPH_ASSERT(bipartite);
    IGRAPH_ASSERT(igraph_vector_bool_size(&types) == igraph_vcount(&graph));
    igraph_destroy(&graph);

    /* Directed path */
    igraph_small(&graph, 0, IGRAPH_DIRECTED,
                 0,1, 1,2,
                 -1);

    igraph_is_bipartite(&graph, &bipartite, &types);
    IGRAPH_ASSERT(bipartite);
    IGRAPH_ASSERT(igraph_vector_bool_size(&types) == igraph_vcount(&graph));

    /* Odd directed cycle */
    igraph_add_edge(&graph, 2, 0);

    igraph_is_bipartite(&graph, &bipartite, &types);
    IGRAPH_ASSERT(! bipartite);
    IGRAPH_ASSERT(igraph_vector_bool_size(&types) == igraph_vcount(&graph));

    igraph_destroy(&graph);

    igraph_vector_bool_destroy(&types);

    VERIFY_FINALLY_STACK();

    return 0;
}
