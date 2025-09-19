
#include <igraph.h>
#include <float.h>

int main(void) {
    igraph_t graph;
    igraph_vector_t pagerank;
    igraph_real_t value;

    /* Initialize the library. */
    igraph_setup();

    /* Create a directed graph */
    igraph_kautz(&graph, 2, 3);

    /* Initialize the vector where the results will be stored */
    igraph_vector_init(&pagerank, 0);

    igraph_pagerank(&graph, /* weights */ NULL,
                    &pagerank, &value,
            /* damping */ 0.85, IGRAPH_DIRECTED,
                    igraph_vss_all(), IGRAPH_PAGERANK_ALGO_PRPACK,
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
