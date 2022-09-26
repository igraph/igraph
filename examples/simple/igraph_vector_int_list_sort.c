
#include <igraph.h>
#include <stdio.h>

int main(void) {
    igraph_t graph;
    igraph_vector_int_list_t cliques;
    igraph_integer_t i, n;

    /* Set a random seed to make the program deterministic */
    igraph_rng_seed(igraph_rng_default(), 31415);

    /* Create a random graph with a given number of vertices and edges */
    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 15, 80, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

    /* Find all maximal cliques in the graph */
    igraph_vector_int_list_init(&cliques, 0);
    igraph_maximal_cliques(&graph, &cliques, -1, -1);

    /* Print the cliques in lexicographical order */
    printf("Maximal cliques in lexicographical order:\n");
    igraph_vector_int_list_sort(&cliques, igraph_vector_int_lex_cmp);
    n = igraph_vector_int_list_size(&cliques);
    for (i=0; i < n; ++i) {
        igraph_vector_int_print(igraph_vector_int_list_get_ptr(&cliques, i));
    }

    /* Print the cliques in colexicographical order */
    printf("\nMaximal cliques in colexicographical order:\n");
    igraph_vector_int_list_sort(&cliques, igraph_vector_int_colex_cmp);
    n = igraph_vector_int_list_size(&cliques);
    for (i=0; i < n; ++i) {
        igraph_vector_int_print(igraph_vector_int_list_get_ptr(&cliques, i));
    }

    /* Destroy data structures when we no longer need them */
    igraph_vector_int_list_destroy(&cliques);
    igraph_destroy(&graph);

    return 0;
}
