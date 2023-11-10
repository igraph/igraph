#include <igraph.h>
#include "test_utilities.h"


int main(void) {

    igraph_t g;

    printf("The graph from for dispersion test:\n");
    igraph_small(&g, 12, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2,
                 1, 2, 1, 3, 1, 4, 1, 5,
                 2, 3, 2, 5, 2, 7,
                 3, 5,
                 4, 5,
                 5, 7,
                 7, 9, 7, 10,
                 8, 9, 8, 10,
                 9, 10,
                 11, 0, 11, 1, 11, 2, 11, 3, 11, 4, 11, 5, 11, 6, 11, 7, 11, 8, 11, 9, 11, 10,
                 -1);

    IGRAPH_ASSERT(igraph_dispersion(&g, 11, 7) == 4);
    IGRAPH_ASSERT(igraph_dispersion(&g, 11, 1) == 1);

    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
