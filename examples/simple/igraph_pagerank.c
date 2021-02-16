
#include <igraph.h>
#include <float.h>

int main() {
    igraph_t graph;
    igraph_vector_t pagerank;
    igraph_real_t value;

    /* Create a directed graph */
    igraph_kautz(&graph, 2, 3);

    /* Initialize the vector where the results will be stored */
    igraph_vector_init(&pagerank, 0);

    igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_PRPACK,
                    &pagerank, &value,
                    igraph_vss_all(), IGRAPH_DIRECTED,
                    /* damping */ 0.85, /* weights */ NULL,
                    NULL /* not needed with PRPACK method */);

    /* Check that the eigenvalue is 1, as expected. */
    if (fabs(value - 1.0) > 32*DBL_EPSILON) {
        fprintf(stderr, "PageRank failed to converge.\n");
        return 1;
    }

    /* Output the result */
    igraph_vector_print(&pagerank);

    /* Destroy data structure when no longer needed */
    igraph_vector_destroy(&pagerank);
    igraph_destroy(&graph);

    return 0;
}
