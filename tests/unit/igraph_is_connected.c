
#include <igraph.h>

#include "test_utilities.h"

int main() {
    igraph_t graph;
    igraph_bool_t conn;

    /* Null graph */
    igraph_empty(&graph, 0, IGRAPH_DIRECTED);

    igraph_is_connected(&graph, &conn, IGRAPH_WEAK);
    IGRAPH_ASSERT(! conn);

    igraph_is_connected(&graph, &conn, IGRAPH_STRONG);
    IGRAPH_ASSERT(! conn);

    igraph_destroy(&graph);

    /* Singleton graph */
    igraph_empty(&graph, 1, IGRAPH_DIRECTED);

    igraph_is_connected(&graph, &conn, IGRAPH_WEAK);
    IGRAPH_ASSERT(conn);

    igraph_is_connected(&graph, &conn, IGRAPH_STRONG);
    IGRAPH_ASSERT(conn);

    igraph_destroy(&graph);

    /* Two isolated vertices, one with a self-loop */
    igraph_small(&graph, 2, IGRAPH_DIRECTED,
                 0,0, -1);

    igraph_is_connected(&graph, &conn, IGRAPH_WEAK);
    IGRAPH_ASSERT(! conn);

    igraph_is_connected(&graph, &conn, IGRAPH_STRONG);
    IGRAPH_ASSERT(! conn);

    igraph_destroy(&graph);

    /* Two isolated vertices, three self-loops */
    igraph_small(&graph, 2, IGRAPH_DIRECTED,
                 0,0, 0,0, 1,1, -1);

    igraph_is_connected(&graph, &conn, IGRAPH_WEAK);
    IGRAPH_ASSERT(! conn);

    igraph_is_connected(&graph, &conn, IGRAPH_STRONG);
    IGRAPH_ASSERT(! conn);

    igraph_destroy(&graph);

    /* Weakly connected directed */
    igraph_small(&graph, 4, IGRAPH_DIRECTED,
                 0,1, 2,0, 1,2, 3,2,
                 -1);

    igraph_is_connected(&graph, &conn, IGRAPH_WEAK);
    IGRAPH_ASSERT(conn);

    igraph_is_connected(&graph, &conn, IGRAPH_STRONG);
    IGRAPH_ASSERT(! conn);

    igraph_destroy(&graph);

    /* Directed cycle */
    igraph_small(&graph, 4, IGRAPH_DIRECTED,
                 0,1, 2,0, 1,3, 3,2,
                 -1);

    igraph_is_connected(&graph, &conn, IGRAPH_WEAK);
    IGRAPH_ASSERT(conn);

    igraph_is_connected(&graph, &conn, IGRAPH_STRONG);
    IGRAPH_ASSERT(conn);

    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
