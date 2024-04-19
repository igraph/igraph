
#include <igraph.h>

#include "test_utilities.h"

int main(void) {

    igraph_t graph;
    igraph_bool_t bipartite, acyclic, has_loop;
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

    /* Singleton with self-loop */
    igraph_small(&graph, 1, IGRAPH_UNDIRECTED,
                 0, 0,
                 -1);
    igraph_is_bipartite(&graph, &bipartite, &types);
    IGRAPH_ASSERT(! bipartite);
    IGRAPH_ASSERT(igraph_vector_bool_size(&types) == igraph_vcount(&graph));

    /* Test cache usage */
    igraph_has_loop(&graph, &has_loop);
    igraph_is_bipartite(&graph, &bipartite, NULL);
    IGRAPH_ASSERT(! bipartite);

    igraph_destroy(&graph);

    /* Directed path */
    igraph_small(&graph, 0, IGRAPH_DIRECTED,
                 0,1, 1,2,
                 -1);

    igraph_is_bipartite(&graph, &bipartite, &types);
    IGRAPH_ASSERT(bipartite);
    IGRAPH_ASSERT(igraph_vector_bool_size(&types) == igraph_vcount(&graph));

    /* Test cache usage */
    bipartite = false;
    igraph_is_forest(&graph, &acyclic, NULL, IGRAPH_ALL);
    igraph_is_bipartite(&graph, &bipartite, NULL);
    IGRAPH_ASSERT(bipartite);

    /* Odd directed cycle */
    igraph_add_edge(&graph, 2, 0);
    igraph_invalidate_cache(&graph);

    igraph_is_bipartite(&graph, &bipartite, &types);
    IGRAPH_ASSERT(! bipartite);
    IGRAPH_ASSERT(igraph_vector_bool_size(&types) == igraph_vcount(&graph));

    /* Test cache usage */
    bipartite = true;
    igraph_is_forest(&graph, &acyclic, NULL, IGRAPH_ALL);
    igraph_is_bipartite(&graph, &bipartite, NULL);
    IGRAPH_ASSERT(! bipartite);

    igraph_destroy(&graph);

    igraph_vector_bool_destroy(&types);

    VERIFY_FINALLY_STACK();

    return 0;
}
