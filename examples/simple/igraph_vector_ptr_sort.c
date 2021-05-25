
#include <igraph.h>
#include <stdio.h>

int main() {
    igraph_t graph;
    igraph_vector_ptr_t cliques;
    long int i, n;

    /* Set a random seed to make the program deterministic */
    igraph_rng_seed(igraph_rng_default(), 31415);

    /* Create a random graph with a given number of vertices and edges */
    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 15, 80, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

    /* Find all maximal cliques in the graph */
    igraph_vector_ptr_init(&cliques, 0);
    igraph_maximal_cliques(&graph, &cliques, -1, -1);

    /* Print the cliques in lexicographical order */
    printf("Maximal cliques in lexicographical order:\n");
    igraph_vector_ptr_sort(&cliques, igraph_vector_lex_cmp);
    n = igraph_vector_ptr_size(&cliques);
    for (i=0; i < n; ++i) {
        igraph_vector_print(VECTOR(cliques)[i]);
    }

    /* Print the cliques in colexicographical order */
    printf("\nMaximal cliques in colexicographical order:\n");
    igraph_vector_ptr_sort(&cliques, igraph_vector_colex_cmp);
    n = igraph_vector_ptr_size(&cliques);
    for (i=0; i < n; ++i) {
        igraph_vector_print(VECTOR(cliques)[i]);
    }

    /* Destroy data structures when we no longer need them */

    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&cliques, igraph_vector_destroy);
    igraph_vector_ptr_destroy_all(&cliques);

    igraph_destroy(&graph);

    return 0;
}
