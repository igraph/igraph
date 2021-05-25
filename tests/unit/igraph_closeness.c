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

    igraph_closeness(&g, &vector_actual_results /*store results here*/,
                      NULL, NULL,
                      /*calculating for all vectors in the graph*/ igraph_vss_all(),
                      IGRAPH_ALL  /*graph is "undirected"*/,
                      NULL /*unweighted*/, /*not normalised*/ 0);

    printf("Non normalised results below\n");
    print_vector(&vector_actual_results);

    /* NORMALISED TEST BELOW */

    igraph_closeness(&g, &vector_actual_results /*store results here*/,
                      NULL, NULL,
                      /*calculating for all vectors in the graph*/ igraph_vss_all(),
                      IGRAPH_ALL  /*graph is "undirected"*/,
                      NULL, /*normalised*/ 1);

    printf("\nNormalised results below\n");
    print_vector(&vector_actual_results);

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

    igraph_closeness(&g, &vector_actual_results /*store results here*/,
                      NULL, NULL,
                      /*calculating for all vectors in the graph*/ igraph_vss_all(),
                      IGRAPH_ALL  /*graph is "undirected"*/,
                      &vector_weights, /*not normalised*/ 0);

    printf("Non normalised test below\n");

    print_vector(&vector_actual_results);

    /* NORMALISED TEST BELOW */

    printf("\nNormalised test below\n");

    igraph_closeness(&g, &vector_actual_results /*store results here*/,
                      NULL, NULL,
                      /*calculating for all vectors in the graph*/ igraph_vss_all(),
                      IGRAPH_ALL  /*graph is "undirected"*/,
                      &vector_weights, /*normalised*/ 1);

    print_vector(&vector_actual_results);

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

    igraph_closeness(&g, &vector_actual_results /*store results here*/,
                      NULL, NULL,
                      /*calculating for all vectors in the graph*/ igraph_vss_all(),
                      IGRAPH_ALL  /*graph is "undirected"*/,
                      NULL, /*not normalised*/ 0);

    print_vector(&vector_actual_results);

    /* NORMALISED TEST BELOW*/

    printf("\nNormalised test below\n");

    igraph_closeness(&g, &vector_actual_results /*store results here*/,
                      NULL, NULL,
                      /*calculating for all vectors in the graph*/ igraph_vss_all(),
                      IGRAPH_ALL  /*graph is "undirected"*/,
                      NULL, /*normalised*/ 1);

    print_vector(&vector_actual_results);

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

    igraph_closeness(&g, &vector_actual_results /*store results here*/,
                      NULL, NULL,
                      /*calculating for all vectors in the graph*/ igraph_vss_all(),
                      IGRAPH_ALL  /*graph is "undirected"*/,
                      &vector_weights, /*not normalised*/ 0);

    print_vector(&vector_actual_results);

    /* TEST FOR DIRECTED GRAPH
    OUT means the min distance from the curr node to the other node */

    printf("\nDirected graph test below for OUT\n");

    igraph_closeness(&g, &vector_actual_results /*store results here*/,
                      NULL, NULL,
                      /*calculating for all vectors in the graph*/ igraph_vss_all(),
                      IGRAPH_OUT  /*graph is "out directed"*/,
                      &vector_weights, /*not normalised*/ 0);

    print_vector(&vector_actual_results);

    /* IN means the min distance from a node to the curr node */

    printf("\nDirected graph test below for IN\n");

    igraph_closeness(&g, &vector_actual_results /*store results here*/,
                      NULL, NULL,
                      /*calculating for all vectors in the graph*/ igraph_vss_all(),
                      IGRAPH_IN  /*graph is "in directed"*/,
                      &vector_weights, /*not normalised*/ 0);

    print_vector(&vector_actual_results);

    igraph_vector_destroy(&vector_actual_results);
    igraph_destroy(&g);
}

void test_cutoff() {

    igraph_t g;
    igraph_vector_t closeness, reachable;
    igraph_bool_t all_reachable;
    size_t i;
    igraph_real_t cutoff_vec[] = { -1.0, 0.0, 1.0, 2.9, 3.0, 3.1 };

    printf("\n\nUnweighted undirected with cutoff\n");

    igraph_ring(&g, 4, IGRAPH_UNDIRECTED, 0, 0);

    igraph_vector_init(&closeness, 0);
    igraph_vector_init(&reachable, 0);

    for (i=0; i < sizeof(cutoff_vec) / sizeof(igraph_real_t); ++i) {
        printf("\nRange-limited closeness with cutoff %g\n", cutoff_vec[i]);
        igraph_closeness_cutoff(&g, &closeness, &reachable, &all_reachable,
                                igraph_vss_all(), IGRAPH_ALL, NULL, /* normalized */ 1,
                                cutoff_vec[i]);
        printf("Closeness: ");
        print_vector(&closeness);
        printf("Reachable: ");
        print_vector_round(&reachable);
        printf("All reachable: %s\n", all_reachable ? "true" : "false");
    }

    igraph_vector_destroy(&reachable);
    igraph_vector_destroy(&closeness);

    igraph_destroy(&g);
}

void test_cutoff_directed() {

    igraph_t g;
    igraph_vector_t closeness, reachable;
    igraph_bool_t all_reachable;
    size_t i;
    igraph_real_t cutoff_vec[] = { -1.0, 0.0, 1.0, 2.9, 3.0, 3.1 };

    printf("\n\nUnweighted directed with cutoff\n");

    igraph_ring(&g, 4, IGRAPH_DIRECTED, 0, 0);

    igraph_vector_init(&closeness, 0);
    igraph_vector_init(&reachable, 0);

    for (i=0; i < sizeof(cutoff_vec) / sizeof(igraph_real_t); ++i) {
        printf("\nRange-limited directed closeness with cutoff %g\n", cutoff_vec[i]);
        igraph_closeness_cutoff(&g, &closeness, &reachable, &all_reachable,
                                igraph_vss_all(), IGRAPH_OUT, NULL, /* normalized */ 1,
                                cutoff_vec[i]);
        printf("Closeness: ");
        print_vector(&closeness);
        printf("Reachable: ");
        print_vector_round(&reachable);
        printf("All reachable: %s\n", all_reachable ? "true" : "false");
    }

    igraph_vector_destroy(&reachable);
    igraph_vector_destroy(&closeness);

    igraph_destroy(&g);
}

void test_cutoff_weighted() {

    igraph_t g;
    igraph_vector_t closeness, reachable;
    igraph_bool_t all_reachable;
    size_t i;
    igraph_real_t cutoff_vec[] = { -1.0, 0.0, 1.0, 2.9, 3.0, 5.0, 6.0 };
    igraph_vector_t weights;

    printf("\n\nWeighted undirected with cutoff\n");

    igraph_ring(&g, 4, IGRAPH_UNDIRECTED, 0, 0);

    igraph_vector_init(&closeness, 0);
    igraph_vector_init(&reachable, 0);
    igraph_vector_init_seq(&weights, 1, 3);

    for (i=0; i < sizeof(cutoff_vec) / sizeof(igraph_real_t); ++i) {
        printf("\nRange-limited weighted closeness with cutoff %g\n", cutoff_vec[i]);
        igraph_closeness_cutoff(&g, &closeness, &reachable, &all_reachable,
                                igraph_vss_all(), IGRAPH_ALL, &weights, /* normalized */ 1,
                                cutoff_vec[i]);
        printf("Closeness: ");
        print_vector(&closeness);
        printf("Reachable: ");
        print_vector_round(&reachable);
        printf("All reachable: %s\n", all_reachable ? "true" : "false");
    }

    igraph_vector_destroy(&weights);
    igraph_vector_destroy(&reachable);
    igraph_vector_destroy(&closeness);

    igraph_destroy(&g);
}

void test_edge_cases() {

    igraph_t g;
    igraph_vector_t closeness, reachable;
    igraph_bool_t all_reachable;
    int n;

    printf("\n\nEdgeless graphs\n");

    igraph_vector_init(&closeness, 0);
    igraph_vector_init(&reachable, 0);

    for (n=0; n <= 2; ++n) {
        printf("\nEdgeless graph with %d vertices\n", n);
        igraph_empty(&g, n, IGRAPH_UNDIRECTED);

        igraph_closeness(&g, &closeness, &reachable, &all_reachable, igraph_vss_all(), IGRAPH_ALL, NULL, 1);
        printf("Closeness: ");
        print_vector(&closeness);
        printf("Reachable: ");
        print_vector_round(&reachable);
        printf("All reachable: %s\n", all_reachable ? "true" : "false");

        igraph_destroy(&g);
    }

    igraph_vector_destroy(&reachable);
    igraph_vector_destroy(&closeness);

}

int main() {

    simple_test_case_no_weights_undirected();
    simple_test_case_with_weights_undirected();
    advanced_test_case_no_weights_undirected();
    advanced_test_case_with_weights();

    test_cutoff();
    test_cutoff_directed();
    test_cutoff_weighted();

    test_edge_cases();

    VERIFY_FINALLY_STACK();

    return 0;
}
