
#include <igraph.h>

int main(void) {
    igraph_t graph;
    igraph_vector_int_t colors;
    igraph_bool_t valid_coloring;

    /* Setting a seed makes the result of erdos_renyi_game_gnm deterministic. */
    igraph_rng_seed(igraph_rng_default(), 42);

    /* IGRAPH_UNDIRECTED and IGRAPH_NO_LOOPS are both equivalent to 0/FALSE, but
       communicate intent better in this context. */
    igraph_erdos_renyi_game_gnm(&graph, 1000, 10000, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

    /* As with all igraph functions, the vector in which the result is returned must
       be initialized in advance. */
    igraph_vector_int_init(&colors, 0);
    igraph_vertex_coloring_greedy(&graph, &colors, IGRAPH_COLORING_GREEDY_COLORED_NEIGHBORS);


    /* Verify that the colouring is valid, i.e. no two adjacent vertices have the same colour. */
    igraph_is_vertex_coloring(&graph, &colors, &valid_coloring);

    if (! valid_coloring) {
        printf("Inconsistent coloring!\n");
        abort();
    }

    /* Destroy data structure when we are done. */
    igraph_vector_int_destroy(&colors);
    igraph_destroy(&graph);

    return 0;
}
