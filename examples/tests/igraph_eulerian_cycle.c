
#include <stdio.h>
#include <igraph.h>
#include "test_utilities.inc"

int main() {

    igraph_t graph;
    /*igraph_vector_t res;*/
    igraph_vector_t edge_res, vector_res;
    igraph_es_t es;
    igraph_vs_t vs;
/*

    igraph_vector_init(&res, 0);
    igraph_eulerian_cycle(&graph, &res, NULL);
    print_vector_round(&res, stdout);
    igraph_vector_destroy(&res);
    igraph_vector_init(&res, 0);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0, 1, -1);
    igraph_es_1(&es, 0);
    igraph_delete_edges(&graph, es);
    igraph_vs_1(&vs, 1);
    igraph_delete_vertices(&graph, vs);
    igraph_vs_1(&vs, 0);
    igraph_delete_vertices(&graph, vs);
    igraph_eulerian_cycle(&graph, &res, NULL);
    print_vector_round(&res, stdout);

    igraph_vector_destroy(&res);
    igraph_vector_init(&res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,3, 3,4, 4,5, 5,2,
                2,6, 6,4, 4,8, 2,8, 2,7, 0,7, -1);
    igraph_eulerian_cycle(&graph, &res, NULL);
    print_vector_round(&res, stdout);

    igraph_vector_destroy(&res);
    igraph_vector_init(&res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,0, 0,0, 0,0, -1);
    igraph_eulerian_cycle(&graph, &res, NULL);
    print_vector_round(&res, stdout);

    igraph_vector_destroy(&res);
    igraph_vector_init(&res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0, 1, -1);
    igraph_es_1(&es, 0);
    igraph_delete_edges(&graph, es);
    igraph_vs_1(&vs, 1);
    igraph_delete_vertices(&graph, vs);
    igraph_eulerian_cycle(&graph, &res, NULL);
    print_vector_round(&res, stdout);

    igraph_destroy(&graph);
    igraph_vector_destroy(&res);
    igraph_vector_init(&res, 0);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0, 1, -1);
    igraph_es_1(&es, 0);
    igraph_delete_edges(&graph, es);
    igraph_vs_1(&vs, 1);
    igraph_delete_vertices(&graph, vs);
    igraph_vs_1(&vs, 0);
    igraph_delete_vertices(&graph, vs);
    igraph_eulerian_cycle(&graph, &res, NULL);
    print_vector_round(&res, stdout);

    igraph_vector_destroy(&res);
    igraph_vector_init(&res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1, -1);
    igraph_es_1(&es, 0);
    igraph_delete_edges(&graph, es);
    igraph_eulerian_cycle(&graph, &res, NULL);
    print_vector_round(&res, stdout);

    igraph_destroy(&graph);
    igraph_vector_destroy(&res);
    igraph_vector_init(&res, 0);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1, 1,2, 2,0, -1);
    igraph_eulerian_cycle(&graph, &res, NULL);
    print_vector_round(&res, stdout);

    igraph_vector_destroy(&res);
    igraph_vector_init(&res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,3, 3,4, 4,0, 0,2, 2,1, 1,0, -1);
    igraph_eulerian_cycle(&graph, &res, NULL);
    print_vector_round(&res, stdout);

    igraph_vector_destroy(&res);
    igraph_vector_init(&res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,6, 6,4, 4,5, 5,0, 0,1, 1,2,
                2,3, 3,4, 4,2, 2,0, -1);
    igraph_eulerian_cycle(&graph, &res, NULL);
    print_vector_round(&res, stdout);

    igraph_vector_destroy(&res);
    igraph_vector_init(&res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0, 1, -1);
    igraph_es_1(&es, 0);
    igraph_delete_edges(&graph, es);
    igraph_eulerian_cycle(&graph, &res, NULL);
    print_vector_round(&res, stdout);

    igraph_destroy(&graph);
    igraph_vector_destroy(&res);
    igraph_vector_init(&res, 0);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0, 1, -1);
    igraph_es_1(&es, 0);
    igraph_delete_edges(&graph, es);
    igraph_vs_1(&vs, 1);
    igraph_delete_vertices(&graph, vs);
    igraph_eulerian_cycle(&graph, &res, NULL);
    print_vector_round(&res, stdout);

    igraph_destroy(&graph);
    igraph_vector_destroy(&res);

*/

    //printf("hello\n");

    igraph_vector_init(&edge_res, 0);
    igraph_vector_init(&vector_res, 0);
    //printf("hello 1.5\n");
    igraph_eulerian_cycle(&graph, &edge_res, &vector_res);
    //printf("hello 1.6\n");
    print_vector_round(&edge_res, stdout);
    print_vector_round(&vector_res, stdout);
    printf("\n");

    //printf("hello 2\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0, 1, -1);
    igraph_es_1(&es, 0);
    igraph_delete_edges(&graph, es);
    igraph_vs_1(&vs, 1);
    igraph_delete_vertices(&graph, vs);
    igraph_vs_1(&vs, 0);
    igraph_delete_vertices(&graph, vs);
    igraph_eulerian_cycle(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res, stdout);
    print_vector_round(&edge_res, stdout);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,3, 3,4, 4,5, 5,2,
                2,6, 6,4, 4,8, 2,8, 2,7, 0,7, -1);
    igraph_eulerian_cycle(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res, stdout);
    print_vector_round(&vector_res, stdout);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,0, 0,0, 0,0, -1);
    igraph_eulerian_cycle(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res, stdout);
    print_vector_round(&vector_res, stdout);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0, 1, -1);
    igraph_es_1(&es, 0);
    igraph_delete_edges(&graph, es);
    igraph_vs_1(&vs, 1);
    igraph_delete_vertices(&graph, vs);
    igraph_eulerian_cycle(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res, stdout);
    print_vector_round(&vector_res, stdout);
    printf("\n");

    igraph_destroy(&graph);
    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0, 1, -1);
    igraph_es_1(&es, 0);
    igraph_delete_edges(&graph, es);
    igraph_vs_1(&vs, 1);
    igraph_delete_vertices(&graph, vs);
    igraph_vs_1(&vs, 0);
    igraph_delete_vertices(&graph, vs);
    igraph_eulerian_cycle(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res, stdout);
    print_vector_round(&vector_res, stdout);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1, -1);
    igraph_es_1(&es, 0);
    igraph_delete_edges(&graph, es);
    igraph_eulerian_cycle(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res, stdout);
    print_vector_round(&vector_res, stdout);
    printf("\n");

    igraph_destroy(&graph);
    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1, 1,2, 2,0, -1);
    igraph_eulerian_cycle(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res, stdout);
    print_vector_round(&vector_res, stdout);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,3, 3,4, 4,0, 0,2, 2,1, 1,0, -1);
    igraph_eulerian_cycle(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res, stdout);
    print_vector_round(&vector_res, stdout);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,6, 6,4, 4,5, 5,0, 0,1, 1,2,
                2,3, 3,4, 4,2, 2,0, -1);
    igraph_eulerian_cycle(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res, stdout);
    print_vector_round(&vector_res, stdout);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0, 1, -1);
    igraph_es_1(&es, 0);
    igraph_delete_edges(&graph, es);
    igraph_eulerian_cycle(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res, stdout);
    print_vector_round(&vector_res, stdout);
    printf("\n");

    igraph_destroy(&graph);
    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vector_res);
    igraph_vector_init(&vector_res, 0);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0, 1, -1);
    igraph_es_1(&es, 0);
    igraph_delete_edges(&graph, es);
    igraph_vs_1(&vs, 1);
    igraph_delete_vertices(&graph, vs);
    igraph_eulerian_cycle(&graph, &edge_res, &vector_res);
    print_vector_round(&edge_res, stdout);
    print_vector_round(&vector_res, stdout);

    igraph_destroy(&graph);
    igraph_vector_destroy(&edge_res);
    igraph_vector_destroy(&vector_res);

    VERIFY_FINALLY_STACK();
    
    return 0;
}
