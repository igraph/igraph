#include <stdio.h>
#include <igraph.h>
#include "test_utilities.inc"

int main() {

    igraph_t graph;
    igraph_vector_t res;

    igraph_vector_init(&res, 0);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,3, 3,4, 4,5, 5,2,
                2,6, 6,4, 4,8, 2,8, 2,7, 0,7, -1);
    igraph_eulerian_cycle(&graph, &res);
    print_vector_round(&res, stdout);

    igraph_vector_destroy(&res);
    igraph_vector_init(&res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1, 1,2, 2,0, -1);
    igraph_eulerian_cycle(&graph, &res);
    print_vector_round(&res, stdout);

    igraph_vector_destroy(&res);
    igraph_vector_init(&res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,3, 3,4, 4,0, 0,2, 2,1, 1,0, -1);
    igraph_eulerian_cycle(&graph, &res);
    print_vector_round(&res, stdout);

    igraph_vector_destroy(&res);
    igraph_vector_init(&res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,6, 6,4, 4,5, 5,0, 0,1, 1,2,
                2,3, 3,4, 4,2, 2,0, -1);
    igraph_eulerian_cycle(&graph, &res);
    print_vector_round(&res, stdout);

    igraph_destroy(&graph);
    igraph_vector_destroy(&res);
    
    return 0;
}