
#include <igraph.h>
#include <stdio.h>
#include "test_utilities.inc"

int main() {

    igraph_t graph;
    igraph_bool_t has_path, has_cycle;

    /* undirected cases */
    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,1 , 1,2, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 1,2 , 2,3, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 1,2 , 2,3, 3,1, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,1, 0,1,0,1,-1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,0, 0,0,0,0,-1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,3, 3,4, 4,5, 5,2,
                2,6, 6,4, 4,8, 2,8, 2,7, 0,7, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,1 , 1,2, 2,3, 2,4 , 3,5 , 4,5,
                4,6, 0,6, 6,7, 1,7, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,1 , 1,2, 2,3, 3,4 , 2,4 , 1,5,
                0,5 , -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,4, 3,4, 1,3, 2,5, 4,5, 2,6, 1,6, 0,4, 6,5, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 1,3 , 3,5 , 5,6 , 6,3 , 3,1 , 1,2 , 
                                        2,2 , 2,4 , 2,4 , 4,3 , 3,2 , 4,6 , -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED,  7,8, 8,9, 9,7, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED,  0,1, 2,3, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED,  0,1, 2,3, 3,1, 4,5, 5,6, 6,4, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    /* two disconnected self loops */
    igraph_small(&graph, 0, IGRAPH_UNDIRECTED,  1,1, 2,2, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    /* one self loop and one disconnected multiedge selfloop */
    igraph_small(&graph, 0, IGRAPH_UNDIRECTED,  1,1, 1,1, 2,2, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph); 

    /* multiple self-loop singletons */
    igraph_small(&graph, 0, IGRAPH_UNDIRECTED,  0,0 , 1,1 , 1,1 , -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    /* no edges, multiple vertices */
    igraph_small(&graph, 4, IGRAPH_UNDIRECTED, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    /* no edges except one self loop, multiple vertices */
    igraph_small(&graph, 4, IGRAPH_UNDIRECTED, 0,0, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    /* directed cases*/
    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1 , 1,2, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1, 1,2, 2,0, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1 , 1,2, 1,3, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1 , 1,3, 3,2, 2,0 , 2,1, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,3, 3,4, 4,0, 0,2, 2,1, 1,0, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,6, 6,4, 4,5, 5,0, 0,1, 1,2,
                2,3, 3,4, 4,2, 2,0, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    /* multiedges */
    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1, 0,1, 1,2, 2,1, 1,3, 3,4, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    /* one self loop */
    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1, 0,0, 1,2, 2,1, 1,3, 3,4, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    /* multiedges and one self loop */
    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,2 , 2,4 , 4,5 , 5,2 , 2,0 , 0,1 , 
                                        1,1 , 1,3 , 1,3 , 3,2 , 2,1 , 3,5 , -1);    
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 1,3 , 3,5 , 5,6 , 6,3 , 3,1 , 1,2 , 
                                        2,2 , 2,4 , 2,4 , 4,3 , 3,2 , 4,6 , -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1, 0,1 , 0,1, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    /* disconnected graphs and self loops, both undirected and directed */
    /* disconnected with singleton vertices */
    igraph_small(&graph, 0, IGRAPH_DIRECTED, 8,9, 9,10, 10,8, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    /* two disconnected self loops, directed */
    igraph_small(&graph, 0, IGRAPH_DIRECTED,  1,1, 2,2, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    /* one self loop and one disconnected multiedge selfloop, directed */
    igraph_small(&graph, 0, IGRAPH_DIRECTED,  1,1, 1,1, 2,2, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    /* no edges, multiple vertices, directed */
    igraph_small(&graph, 4, IGRAPH_DIRECTED, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    /* no edges except one self loop, multiple vertices, directed */
    igraph_small(&graph, 4, IGRAPH_DIRECTED, 0,0, -1);
    igraph_is_eulerian(&graph, &has_path, &has_cycle);
    printf("%d %d\n", has_path, has_cycle);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();
    
    return 0;
}
