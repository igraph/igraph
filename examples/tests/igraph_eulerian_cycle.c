#include <stdio.h>
#include "test_utilities.inc"

int main() {

    igraph_t graph, graph2, graph3, graph4;
    igraph_vector_t res, res2, res3, res4;

    igraph_vector_init(&res, 0);
    igraph_vector_init(&res2, 0);  
    igraph_vector_init(&res3, 0);  
    igraph_vector_init(&res4, 0);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,3, 3,4, 4,5, 5,2,
                2,6, 6,4, 4,8, 2,8, 2,7, 0,7, -1);
    igraph_eulerian_cycle(&graph, &res);
    print_vector_round(&res, stdout);

    igraph_small(&graph2, 0, IGRAPH_DIRECTED, 0,1, 1,2, 2,0, -1);
    igraph_eulerian_cycle(&graph2, &res2);
    print_vector_round(&res2, stdout);

    igraph_small(&graph3, 0, IGRAPH_DIRECTED, 0,3, 3,4, 4,0, 0,2, 2,1, 1,0, -1);
    igraph_eulerian_cycle(&graph3, &res3);
    print_vector_round(&res3, stdout);

    igraph_small(&graph4, 0, IGRAPH_DIRECTED, 0,6, 6,4, 4,5, 5,0, 0,1, 1,2,
                2,3, 3,4, 4,2, 2,0, -1);
    igraph_eulerian_cycle(&graph4, &res4);
    print_vector_round(&res4, stdout);

    igraph_destroy(&graph);
    igraph_destroy(&graph2);
    igraph_destroy(&graph3);
    igraph_destroy(&graph4);

    igraph_vector_destroy(&res);
    igraph_vector_destroy(&res2);
    igraph_vector_destroy(&res3);
    igraph_vector_destroy(&res4);
    
    return 0;
}