#include <igraph.h>
#include <stdio.h>
#include "test_utilities.inc"


void simple_test_case_no_weights_undirected() {

    igraph_t g;
    igraph_vector_t vector_actual_results;

    printf("Simple test case, no weights, undirected\n");

    igraph_vector_init(&vector_actual_results, 0);

    igraph_small(&g, 0, IGRAPH_DIRECTED, 0,1 , 1,2, -1);

    /* NOT NORMALISED TEST BELOW */

    igraph_closeness_estimate(&g, &vector_actual_results /*store results here*/, 
                             /*calculating for all vectors in the graph*/ igraph_vss_all(), 
                             IGRAPH_ALL  /*graph is "undirected"*/, /*cutoff*/ 0 /*calculating exact centrality*/, 
                             NULL /*unweighted*/, /*not normalised*/ 0);

    printf("Non normalised results below\n");
    print_vector(&vector_actual_results, stdout);

    /* NORMALISED TEST BELOW */

    igraph_closeness_estimate(&g, &vector_actual_results /*store results here*/, 
                             /*calculating for all vectors in the graph*/ igraph_vss_all(), 
                             IGRAPH_ALL  /*graph is "undirected"*/, /*cutoff*/ 0 /*calculating exact centrality*/, 
                             NULL, /*normalised*/ 1);

    printf("\nNormalised results below\n");
    print_vector(&vector_actual_results, stdout); 

    igraph_vector_destroy(&vector_actual_results);
    igraph_destroy(&g);
}

void simple_test_case_with_weights_undirected() {

    igraph_t g;
    igraph_vector_t vector_edges, vector_weights, vector_actual_results;

    igraph_real_t real_edges[] = {0,1 , 1,2};
    igraph_real_t real_weights[] = {3, 5};

    printf("\nSimple test case, with weights, undirected\n");

    igraph_vector_init(&vector_actual_results, 0);
    igraph_vector_view(&vector_edges, real_edges, sizeof(real_edges)/sizeof(igraph_real_t));
    igraph_create(&g, &vector_edges, /*number of vertices*/ 2, IGRAPH_DIRECTED);
    
    igraph_vector_view(&vector_weights, real_weights, 
                       sizeof(real_weights)/sizeof(igraph_real_t));

    /* NOT NORMALISED TEST BELOW */

    igraph_closeness_estimate(&g, &vector_actual_results /*store results here*/, 
                             /*calculating for all vectors in the graph*/ igraph_vss_all(), 
                             IGRAPH_ALL  /*graph is "undirected"*/, /*cutoff*/ 0 /*calculating exact centrality*/, 
                             &vector_weights, /*not normalised*/ 0);
 
    printf("Non normalised test below\n"); 
    
    print_vector(&vector_actual_results, stdout);

    /* NORMALISED TEST BELOW */

    printf("\nNormalised test below\n");

    igraph_closeness_estimate(&g, &vector_actual_results /*store results here*/, 
                             /*calculating for all vectors in the graph*/ igraph_vss_all(), 
                             IGRAPH_ALL  /*graph is "undirected"*/, /*cutoff*/ 0 /*calculating exact centrality*/, 
                             &vector_weights, /*normalised*/ 1);

    print_vector(&vector_actual_results, stdout);

    igraph_vector_destroy(&vector_actual_results);
    igraph_destroy(&g);
}

void advanced_test_case_no_weights_undirected() {

    igraph_t g;
    igraph_vector_t vector_actual_results;

    /* note, denominatory calculated as (shortest dist)*(n-1) for no normalisation
    normalisation excludes n-1 as part of the denominator */

    printf("\nAdvanced test case, no weights, undirected\n");

    igraph_vector_init(&vector_actual_results, 0);   

    igraph_small(&g, 0, IGRAPH_DIRECTED, 1,0 , 0,5 , 5,6 , 5,4, 4,1 , 1,2 , 2,4 , 4,6 ,
                                 2,3 , 3,7 , 7,6 , 2,6 , -1);

    /* NOT NORMALISED TEST BELOW*/

    printf("Non normalised test below\n");

    igraph_closeness_estimate(&g, &vector_actual_results /*store results here*/, 
                             /*calculating for all vectors in the graph*/ igraph_vss_all(), 
                             IGRAPH_ALL  /*graph is "undirected"*/, /*cutoff*/ 0 /*calculating exact centrality*/, 
                             NULL, /*not normalised*/ 0);

    print_vector(&vector_actual_results, stdout);
 
    /* NORMALISED TEST BELOW*/

    printf("\nNormalised test below\n");

    igraph_closeness_estimate(&g, &vector_actual_results /*store results here*/, 
                             /*calculating for all vectors in the graph*/ igraph_vss_all(), 
                             IGRAPH_ALL  /*graph is "undirected"*/, /*cutoff*/ 0 /*calculating exact centrality*/, 
                             NULL, /*normalised*/ 1);

    print_vector(&vector_actual_results, stdout);

    igraph_vector_destroy(&vector_actual_results);
    igraph_destroy(&g);
}

void advanced_test_case_with_weights() {

    igraph_t g;
    igraph_vector_t vector_edges, vector_weights, vector_actual_results;

    igraph_real_t real_edges[] = {1,0 , 0,5 , 5,6 , 5,4, 4,1 , 1,2 , 2,4 , 4,6 ,
                                 2,3 , 3,7 , 7,6 , 6,2};
    
    igraph_real_t real_weights[] = {4, 9, 2, 2, 2, 3, 1, 1, 8, 7, 5, 5};

    printf("\nAdvanced test case, with weights\n");


    igraph_vector_init(&vector_actual_results, 0);
    igraph_vector_view(&vector_edges, real_edges, sizeof(real_edges)/sizeof(igraph_real_t));
    igraph_create(&g, &vector_edges, /*number of vertices*/ 2, IGRAPH_DIRECTED);

    igraph_vector_view(&vector_weights, real_weights, sizeof(real_weights)/sizeof(igraph_real_t));

    /* TEST FOR UNDIRECTED GRAPH */

    printf("Undirected graph test below\n");

    igraph_closeness_estimate(&g, &vector_actual_results /*store results here*/, 
                             /*calculating for all vectors in the graph*/ igraph_vss_all(), 
                             IGRAPH_ALL  /*graph is "undirected"*/, /*cutoff*/ 0 /*calculating exact centrality*/, 
                             &vector_weights, /*not normalised*/ 0);

    print_vector(&vector_actual_results, stdout);

    /* TEST FOR DIRECTED GRAPH
    OUT means the min distance from the curr node to the other node */

    printf("\nDirected graph test below for OUT\n");

    igraph_closeness_estimate(&g, &vector_actual_results /*store results here*/, 
                             /*calculating for all vectors in the graph*/ igraph_vss_all(), 
                             IGRAPH_OUT  /*graph is "out directed"*/, /*cutoff*/ 0 /*calculating exact centrality*/, 
                             &vector_weights, /*not normalised*/ 0);

    print_vector(&vector_actual_results, stdout);

    /* IN means the min distance from a node to the curr node */

    printf("\nDirected graph test below for IN\n");

    igraph_closeness_estimate(&g, &vector_actual_results /*store results here*/, 
                             /*calculating for all vectors in the graph*/ igraph_vss_all(), 
                             IGRAPH_IN  /*graph is "in directed"*/, /*cutoff*/ 0 /*calculating exact centrality*/, 
                             &vector_weights, /*not normalised*/ 0);

    print_vector(&vector_actual_results, stdout); 

    igraph_vector_destroy(&vector_actual_results);
    igraph_destroy(&g);
}

int main() {

    simple_test_case_no_weights_undirected();
    simple_test_case_with_weights_undirected();
    advanced_test_case_no_weights_undirected();
    advanced_test_case_with_weights();

    return 0;
}