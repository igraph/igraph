#include <igraph.h>

#include "test_utilities.h"

igraph_error_t verify_coloring(igraph_t *graph, igraph_vector_int_t *colors){
    /* Verify that the colouring is valid, i.e. no two adjacent vertices have the same colour. */
    igraph_integer_t i;
    /* Store the edge count to avoid the overhead from igraph_ecount in the for loop. */
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    for (i = 0; i < no_of_edges; ++i) {
        if ( VECTOR(*colors)[ IGRAPH_FROM(graph, i) ] == VECTOR(*colors)[ IGRAPH_TO(graph, i) ]  ) {
            printf("Inconsistent coloring! Vertices %" IGRAPH_PRId " and %" IGRAPH_PRId " are adjacent but have the same color.\n",
                   IGRAPH_FROM(graph, i), IGRAPH_TO(graph, i));
            return IGRAPH_FAILURE;
        }
    }
    igraph_integer_t vertex_count = igraph_vcount(graph);
    /*verify all assigned colors are valid*/
    for(i = 0 ; i < vertex_count ; i++){
        igraph_integer_t color = VECTOR(*colors)[i];
        if(color<0 || color>=vertex_count){
            printf("The vertex %" IGRAPH_PRId " has invalid color %" IGRAPH_PRId "\n", i, color);
            return IGRAPH_FAILURE;
        }
    }
    return IGRAPH_SUCCESS;
}

int main(void) {
    igraph_t graph;
    igraph_vector_int_t colors;

    /* Set default seed to get reproducible results. */
    igraph_rng_seed(igraph_rng_default(), 42);

    /* Simple undirected graph. */
    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 1000, 10000, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

    /* Initialize the vector into which the function will write its results. */
    igraph_vector_int_init(&colors, 0);
    IGRAPH_CHECK(igraph_vertex_coloring_greedy(&graph, &colors, IGRAPH_COLORING_GREEDY_COLORED_NEIGHBORS));
    IGRAPH_CHECK(verify_coloring(&graph, &colors));

    igraph_vector_int_fill(&colors, 0);
    IGRAPH_CHECK(igraph_vertex_coloring_greedy(&graph, &colors, IGRAPH_COLORING_GREEDY_DSATUR));
    IGRAPH_CHECK(verify_coloring(&graph, &colors));

    /* Destroy data structures when we are done with them. */
    igraph_vector_int_destroy(&colors);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}

