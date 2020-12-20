
#include <stdio.h>
#include <igraph.h>
#include "test_utilities.inc"

int main() {

    igraph_t graph;
    igraph_vector_t edge_res, vertex_res;
    igraph_es_t es;
    igraph_vs_t vs;

    igraph_vector_init(&edge_res, 0);
    igraph_vector_init(&vertex_res, 0);
/*    
    igraph_eulerian_cycle(&graph, &edge_res, &vertex_res);
    print_vector_round(&edge_res);
    print_vector_round(&vertex_res);
*/
    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vertex_res);
    igraph_vector_init(&vertex_res, 0);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0, 1, -1);
    igraph_es_1(&es, 0);
    igraph_delete_edges(&graph, es);
    igraph_vs_1(&vs, 1);
    igraph_delete_vertices(&graph, vs);
    igraph_vs_1(&vs, 0);
    igraph_delete_vertices(&graph, vs);
    igraph_eulerian_cycle(&graph, &edge_res, &vertex_res);
    print_vector_round(&edge_res);
    print_vector_round(&edge_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vertex_res);
    igraph_vector_init(&vertex_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,3, 3,4, 4,5, 5,2,
                2,6, 6,4, 4,8, 2,8, 2,7, 0,7, -1);
    igraph_eulerian_cycle(&graph, &edge_res, &vertex_res);
    print_vector_round(&edge_res);
    print_vector_round(&vertex_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vertex_res);
    igraph_vector_init(&vertex_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,0, 0,0, 0,0, -1);
    igraph_eulerian_cycle(&graph, &edge_res, &vertex_res);
    print_vector_round(&edge_res);
    print_vector_round(&vertex_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vertex_res);
    igraph_vector_init(&vertex_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0, 1, -1);
    igraph_es_1(&es, 0);
    igraph_delete_edges(&graph, es);
    igraph_vs_1(&vs, 1);
    igraph_delete_vertices(&graph, vs);
    igraph_eulerian_cycle(&graph, &edge_res, &vertex_res);
    print_vector_round(&edge_res);
    print_vector_round(&vertex_res);
    printf("\n");

    igraph_destroy(&graph);
    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vertex_res);
    igraph_vector_init(&vertex_res, 0);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0, 1, -1);
    igraph_es_1(&es, 0);
    igraph_delete_edges(&graph, es);
    igraph_vs_1(&vs, 1);
    igraph_delete_vertices(&graph, vs);
    igraph_vs_1(&vs, 0);
    igraph_delete_vertices(&graph, vs);
    igraph_eulerian_cycle(&graph, &edge_res, &vertex_res);
    print_vector_round(&edge_res);
    print_vector_round(&vertex_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vertex_res);
    igraph_vector_init(&vertex_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1, -1);
    igraph_es_1(&es, 0);
    igraph_delete_edges(&graph, es);
    igraph_eulerian_cycle(&graph, &edge_res, &vertex_res);
    print_vector_round(&edge_res);
    print_vector_round(&vertex_res);
    printf("\n");

    igraph_destroy(&graph);
    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vertex_res);
    igraph_vector_init(&vertex_res, 0);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1, 1,2, 2,0, -1);
    igraph_eulerian_cycle(&graph, &edge_res, &vertex_res);
    print_vector_round(&edge_res);
    print_vector_round(&vertex_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vertex_res);
    igraph_vector_init(&vertex_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,3, 3,4, 4,0, 0,2, 2,1, 1,0, -1);
    igraph_eulerian_cycle(&graph, &edge_res, &vertex_res);
    print_vector_round(&edge_res);
    print_vector_round(&vertex_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vertex_res);
    igraph_vector_init(&vertex_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,6, 6,4, 4,5, 5,0, 0,1, 1,2,
                2,3, 3,4, 4,2, 2,0, -1);
    igraph_eulerian_cycle(&graph, &edge_res, &vertex_res);
    print_vector_round(&edge_res);
    print_vector_round(&vertex_res);
    printf("\n");

    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vertex_res);
    igraph_vector_init(&vertex_res, 0);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0, 1, -1);
    igraph_es_1(&es, 0);
    igraph_delete_edges(&graph, es);
    igraph_eulerian_cycle(&graph, &edge_res, &vertex_res);
    print_vector_round(&edge_res);
    print_vector_round(&vertex_res);
    printf("\n");

    igraph_destroy(&graph);
    igraph_vector_destroy(&edge_res);
    igraph_vector_init(&edge_res, 0);
    igraph_vector_destroy(&vertex_res);
    igraph_vector_init(&vertex_res, 0);

    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0, 1, -1);
    igraph_es_1(&es, 0);
    igraph_delete_edges(&graph, es);
    igraph_vs_1(&vs, 1);
    igraph_delete_vertices(&graph, vs);
    igraph_eulerian_cycle(&graph, &edge_res, &vertex_res);
    print_vector_round(&edge_res);
    print_vector_round(&vertex_res);

    igraph_destroy(&graph);
    igraph_vector_destroy(&edge_res);
    igraph_vector_destroy(&vertex_res);

    VERIFY_FINALLY_STACK();
    
    return 0;
}
