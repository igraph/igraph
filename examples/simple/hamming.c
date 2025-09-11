#include <igraph.h>
#include <stdio.h>

int main() {
    igraph_t g;
    igraph_error_t ret;

    // Build Hamming graph H(3,2) → 8 vertices, 12 edges
    ret = igraph_hamming_graph(&g, 3, 2);
    if (ret) {
        fprintf(stderr, "igraph_hamming_graph failed!\n");
        return 1;
    }

    printf("H(3,2):\n");
    printf("  Vertices: %ld\n", (long) igraph_vcount(&g));
    printf("  Edges: %ld\n", (long) igraph_ecount(&g));

    igraph_destroy(&g);
    return 0;
}
