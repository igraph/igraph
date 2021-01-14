
#include <igraph.h>
#include <assert.h>

#include "test_utilities.inc"

int main() {
    igraph_t g1, g2;
    igraph_bool_t res;

    /* undirected graphs */
    igraph_small(&g1, 4, 0,
                 0, 1, 1, 2, 2, 3, 3, 0, -1);
    igraph_small(&g2, 4, 0,
                 1, 0, 1, 2, 2, 3, 3, 0, -1);

    /* a graph is always same as itself */
    res = igraph_is_same_graph(&g1, &g1);
    assert(res);

    /* undirected graphs should be the same no matter
     * the direction of the edges (one is swapped in g2 */
    res = igraph_is_same_graph(&g1, &g2);
    assert(res);

    /* end of undirected */
    igraph_destroy(&g1);
    igraph_destroy(&g2);

    /* directed graphs */
    igraph_small(&g1, 4, 1,
                 0, 1, 1, 2, 2, 3, 3, 0, -1);
    igraph_small(&g2, 4, 1,
                 1, 0, 1, 2, 2, 3, 3, 0, -1);

    /* directed graphs should not be the same if an
     * edge has the opposite direction */
    res = igraph_is_same_graph(&g1, &g2);
    assert(!res);

    igraph_destroy(&g2);

    /* change order of edges, they should be reordered by graph->ii */
    igraph_small(&g2, 4, 1,
                 1, 2, 0, 1, 2, 3, 3, 0, -1);
    res = igraph_is_same_graph(&g1, &g2);
    assert(res);

    /* end of directed */
    igraph_destroy(&g1);
    igraph_destroy(&g2);

    /* undirected vs directed */
    igraph_small(&g1, 4, 0,
                 0, 1, 1, 2, 2, 3, 3, 0, -1);
    igraph_small(&g2, 4, 1,
                 0, 1, 1, 2, 2, 3, 3, 0, -1);
    res = igraph_is_same_graph(&g1, &g2);
    assert(!res);

    /* end of undirected vs directed */
    igraph_destroy(&g1);
    igraph_destroy(&g2);

    VERIFY_FINALLY_STACK();

    return 0;
}
