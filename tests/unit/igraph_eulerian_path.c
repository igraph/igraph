
#include <stdio.h>
#include <igraph.h>
#include "test_utilities.inc"

int main() {

    igraph_t graph;
    igraph_vector_t edge_res, vector_res;
    igraph_es_t es;
    igraph_vs_t vs;

    igraph_vector_init(&edge_res, 0);
    igraph_vector_init(&vector_res, 0);
    /*
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    printf("\n");
    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    */

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0, 1, -1);
    igraph_es_1(&es, 0);
    igraph_delete_edges(&graph, es);
    igraph_vs_1(&vs, 1);
    igraph_delete_vertices(&graph, vs);
    igraph_vs_1(&vs, 0);
    igraph_delete_vertices(&graph, vs);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,1 , 1,2, -1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,0 , -1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,0 ,
                                            1,5, 5,4, 4,1, -1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,1, 1,0, -1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,1 , 1,0 , 1,2, 2,3, 3,4, 4,1 ,-1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,3, 3,4, 4,5, 5,2,
                2,6, 6,4, 4,8, 2,8, 2,7, 0,7, -1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,1 , 1,2, 2,3, 3,4 , 2,4 , 1,5,
                0,5 , -1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,4, 3,4, 1,3, 2,5, 4,5, 2,6, 1,6, 0,4, 6,5, -1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_destroy(&graph);
    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 7,8, 8,9, 9,7, -1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 1,2, 2,3, 3,4, 4,5, 5,6, 6,3,
                3,7, 7,5, 5,9, 3,9, 3,8, 1,8, -1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_destroy(&graph);
    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,1 , 1,2 , 2,0 , 0,3 , 3,2 ,
                                               2,5 , 5,3 , 3,4 , 4,5 , 1,4 , 0,4 ,-1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_destroy(&graph);
    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 2,3 , 3,4 , 4,2 , 2,5 , 5,4 ,
                                               4,7 , 7,5 , 5,6 , 6,7 , 3,6 , 2,6 ,-1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_destroy(&graph);
    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,1 , 1,2, 2,5 , 5,4 , 5,6 , 6,2 ,
                                            2,3 , 3,4 , 4,8, 8,2 , 2,7, 7,0 , -1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_destroy(&graph);
    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED,  0,1, 1,4 , 4,0, 1,2 , 2,3 , 3,5, 5,2, -1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_destroy(&graph);
    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0, 1, -1);
    igraph_es_1(&es, 0);
    igraph_delete_edges(&graph, es);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_destroy(&graph);
    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0, 1, -1);
    igraph_es_1(&es, 0);
    igraph_delete_edges(&graph, es);
    igraph_vs_1(&vs, 1);
    igraph_delete_vertices(&graph, vs);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0, 1, -1);
    igraph_es_1(&es, 0);
    igraph_delete_edges(&graph, es);
    igraph_vs_1(&vs, 1);
    igraph_delete_vertices(&graph, vs);
    igraph_vs_1(&vs, 0);
    igraph_delete_vertices(&graph, vs);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1 , 1,2, -1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,0, 0,0, 0,0, -1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_destroy(&graph);
    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1, 1,2, 2,0, -1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1 , 1,3, 3,2, 2,0 , 2,1, -1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,3, 3,4, 4,0, 0,2, 2,1, 1,0, -1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,6, 6,4, 4,5, 5,0, 0,1, 1,2,
                2,3, 3,4, 4,2, 2,0, -1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1, 3,4, 1,3, 2,1, 1,2, -1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1, 3,4, 1,3, 2,1, 1,2, 0,0, -1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,2 , 2,4 , 4,5 , 5,2 , 2,0 , 0,1 ,
                                        1,1 , 1,3 , 1,3 , 3,2 , 2,1 , 3,5 , -1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 1,3 , 3,5 , 5,6 , 6,3 , 3,1 , 1,2 ,
                                        2,2 , 2,4 , 2,4 , 4,3 , 3,2 , 4,6 , -1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 7,8, 8,9, -1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 7,8, 8,9, 9,7, -1);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0, 1, -1);
    igraph_es_1(&es, 0);
    igraph_delete_edges(&graph, es);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);
    printf("\n");

    igraph_destroy(&graph);
    igraph_vector_destroy(&edge_res);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_vector_init(&edge_res, 0);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0, 1, -1);
    igraph_es_1(&es, 0);
    igraph_delete_edges(&graph, es);
    igraph_vs_1(&vs, 1);
    igraph_delete_vertices(&graph, vs);
    igraph_eulerian_path(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res);
    print_vector_round(&vector_res);

    igraph_destroy(&graph);
    igraph_vector_destroy(&vector_res);
    igraph_vector_destroy(&edge_res);

    VERIFY_FINALLY_STACK();

    return 0;
}
