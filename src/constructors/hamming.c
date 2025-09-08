#include <igraph.h>
#include <math.h>
#include <stdlib.h>

/* 
 * Construct a Hamming graph H(d, q).
 * - d = dimension
 * - q = alphabet size
 * Vertices = q^d
 * Edges connect vertices that differ in exactly 1 position
 */
int igraph_hamming_graph(igraph_t *graph, int d, int q) {
    if (d <= 0 || q <= 1) {
        IGRAPH_ERROR("Both d and q must be positive and q > 1", IGRAPH_EINVAL);
    }

    long n = pow(q, d);  // number of vertices
    igraph_vector_t edges;
    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

    // Represent vertices as numbers 0...(q^d - 1)
    for (long i = 0; i < n; i++) {
        for (long j = i + 1; j < n; j++) {
            int diff = 0;
            long a = i, b = j;

            // Compare base-q digits
            for (int k = 0; k < d; k++) {
                if ((a % q) != (b % q)) diff++;
                a /= q;
                b /= q;
            }

            if (diff == 1) {
                igraph_vector_push_back(&edges, i);
                igraph_vector_push_back(&edges, j);
            }
        }
    }

    // Create graph from edges
    igraph_create(graph, &edges, n, IGRAPH_UNDIRECTED);

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}
