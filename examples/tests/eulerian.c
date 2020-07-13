#include <igraph.h>
#include <stdio.h>

int main(void) {
    igraph_t graph;

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,1 , 1,2, -1);
    assert(is_eulerian(&graph) ==  2);
    print_euler_tour(&graph);

    printf("test 1 done\n");
    
    igraph_t graph2;
    igraph_small(&graph2, 0, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,3, 3,4, 4,5, 5,2,
                2,6, 6,4, 4,8, 2,8, 2,7, 0,7, -1);
    assert(is_eulerian(&graph2) ==  1);
    print_euler_tour(&graph2);
    
    printf("test 2 done\n");

    igraph_t graph3;
    igraph_small(&graph3, 0, IGRAPH_UNDIRECTED, 0,1 , 1,2, 2,3, 2,4 , 3,5 , 4,5,
                4,6, 0,6, 6,7, 1,7, -1);
    assert(is_eulerian(&graph3) ==  0);

    printf("test 3 done\n");

    igraph_t graph4;
    igraph_small(&graph4, 0, IGRAPH_UNDIRECTED, 0,1 , 1,2, 2,3, 3,4 , 2,4 , 1,5,
                0,5 , -1);
    assert(is_eulerian(&graph4) ==  2);
    print_euler_tour(&graph4);

    printf("test 4 done\n");

    igraph_t graph5;
    igraph_small(&graph5, 0, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,4, 3,4, 1,3, 2,5, 4,5, 2,6, 1,6, 0,4, 6,5, -1);
    assert(is_eulerian(&graph5) == 2);
    print_euler_tour(&graph5);

    printf("test 5 done\n");

    printf("all tests done\n");


    return 0;
}