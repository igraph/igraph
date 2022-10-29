
#include <igraph.h>

/* This is a callback function suitable for use with igraph_motifs_randesu_callback().
 * It prints each motif it is calld with. */
igraph_error_t print_motif(const igraph_t *graph, igraph_vector_int_t *vids,
                          igraph_integer_t isoclass, void* extra) {
    printf("Found isoclass %2" IGRAPH_PRId ":  ", isoclass);
    igraph_vector_int_print(vids);
    return IGRAPH_SUCCESS; /* Return 'IGRAPH_SUCCESS': do not interrupt the search. */
}

int main(void) {

    igraph_t graph;
    igraph_vector_t hist;
    igraph_real_t zeros[] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
    igraph_vector_t cut_prob;

    /* Compute the 4-motif distritbuion in Zachary's karate club network. */

    igraph_famous(&graph, "Zachary");
    igraph_vector_init(&hist, 0);

    igraph_motifs_randesu(&graph, &hist, 4, igraph_vector_view(&cut_prob, zeros, 4));

    /* Compute the total number of motifs (connected 4-vertex subgraphs)
     * so that we can print the normalized distribution. */
    igraph_real_t sum = 0.0;
    igraph_integer_t n = igraph_vector_size(&hist);
    for (igraph_integer_t i=0; i < n; i++) {
        if (!isnan(VECTOR(hist)[i])) {
            sum += VECTOR(hist)[i];
        }
    }
    printf("4-motif distribution:\n");
    for (igraph_integer_t i=0; i < n; i++) {
        /* Print NaN values in a platform-independent manner: */
        igraph_real_printf(VECTOR(hist)[i] / sum);
        printf(" ");
    }
    printf("\n\n");

    igraph_vector_destroy(&hist);
    igraph_destroy(&graph);

    /* Identify the vertices of each three-motif in a small Kautz graph. */

    igraph_kautz(&graph, 2, 1);
    igraph_motifs_randesu_callback(&graph, 3, igraph_vector_view(&cut_prob, zeros, 3), &print_motif, NULL);
    igraph_destroy(&graph);

    return 0;
}
