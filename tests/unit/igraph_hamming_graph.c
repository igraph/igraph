#include <igraph.h>
#include "test_utilities.h"

int main(void) {
    igraph_t g;

    IGRAPH_ASSERT(igraph_hamming_graph(&g, 3, 2) == IGRAPH_SUCCESS);

    /* Print something simple and stable to compare in the .out file */
    printf("H(3,2):\n");
    printf("  Vertices: IGRAPH_PRId \n", (igraph_integer_t ) igraph_vcount(&g));
    printf("  Edges: IGRAPH_PRId \n", (igraph_integer_t ) igraph_ecount(&g));

    igraph_destroy(&g);
    VERIFY_FINALLY_STACK();
    return 0;
}
