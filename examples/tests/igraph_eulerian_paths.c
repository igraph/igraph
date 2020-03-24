#include <stdio.h>
#include "test_utilities.inc"

int main() {

    igraph_t graph, graph2, graph3, graph4, graph5, graph6, graph7, graph8, graph9;
    igraph_vector_t res, res2, res3, res4, res5, res6, res7, res8, res9;

    igraph_vector_init(&res, 0);
    igraph_vector_init(&res2, 0);  
    igraph_vector_init(&res3, 0);  
    igraph_vector_init(&res4, 0);
    igraph_vector_init(&res5, 0);
    igraph_vector_init(&res6, 0);  
    igraph_vector_init(&res7, 0);  
    igraph_vector_init(&res8, 0);     
    igraph_vector_init(&res9, 0); 

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,1 , 1,2, -1);
    igraph_eulerian_paths(&graph, &res);
    print_vector(&res, stdout);

    igraph_small(&graph2, 0, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,3, 3,4, 4,5, 5,2,
                2,6, 6,4, 4,8, 2,8, 2,7, 0,7, -1);
    igraph_eulerian_paths(&graph2, &res2);
    print_vector(&res2, stdout);
    
    igraph_small(&graph3, 0, IGRAPH_UNDIRECTED, 0,1 , 1,2, 2,3, 3,4 , 2,4 , 1,5,
                0,5 , -1);
    igraph_eulerian_paths(&graph3, &res3);
    print_vector(&res3, stdout);

    igraph_small(&graph4, 0, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,4, 3,4, 1,3, 2,5, 4,5, 2,6, 1,6, 0,4, 6,5, -1);
    igraph_eulerian_paths(&graph4, &res4);
    print_vector(&res4, stdout);

    igraph_small(&graph5, 0, IGRAPH_DIRECTED, 0,1 , 1,2, -1);
    igraph_eulerian_paths(&graph5, &res5);
    print_vector(&res5, stdout);

    igraph_small(&graph6, 0, IGRAPH_DIRECTED, 0,1, 1,2, 2,0, -1);
    igraph_eulerian_paths(&graph6, &res6);
    print_vector(&res6, stdout);
 
    igraph_small(&graph7, 0, IGRAPH_DIRECTED, 0,1 , 1,3, 3,2, 2,0 , 2,1, -1);
    igraph_eulerian_paths(&graph7, &res7);
    print_vector(&res7, stdout);

    igraph_small(&graph8, 0, IGRAPH_DIRECTED, 0,3, 3,4, 4,0, 0,2, 2,1, 1,0, -1);
    igraph_eulerian_paths(&graph8, &res8);
    print_vector(&res8, stdout);

    igraph_small(&graph9, 0, IGRAPH_DIRECTED, 0,6, 6,4, 4,5, 5,0, 0,1, 1,2,
                2,3, 3,4, 4,2, 2,0, -1);
    igraph_eulerian_paths(&graph9, &res9);
    print_vector(&res9, stdout);

    igraph_destroy(&graph);
    igraph_destroy(&graph2);
    igraph_destroy(&graph3);
    igraph_destroy(&graph4);
    igraph_destroy(&graph5);
    igraph_destroy(&graph6);
    igraph_destroy(&graph7);
    igraph_destroy(&graph8);
    igraph_destroy(&graph9);

    igraph_vector_destroy(&res);
    igraph_vector_destroy(&res2);
    igraph_vector_destroy(&res3);
    igraph_vector_destroy(&res4);
    igraph_vector_destroy(&res5);
    igraph_vector_destroy(&res6);
    igraph_vector_destroy(&res7);
    igraph_vector_destroy(&res8);
    igraph_vector_destroy(&res9);

    return 0;
}