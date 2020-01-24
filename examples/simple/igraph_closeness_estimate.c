#include <igraph.h>
#include <stdio.h>
#include <math.h>

// helper function, not in use anymore but useful for viewing output
void igraph_vector_t_print(igraph_vector_t v, int size) {
    for (int i = 0; i < size; i++) {
        double res = VECTOR(v)[i];
        printf("%lf\n", res);
    }
}

int igraph_check_vectors_equal(igraph_vector_t v, igraph_vector_t w, int size) {
    for (int i = 0; i < size; i++) {
        int v1 = 1000000 * VECTOR(v)[i];
        int v2 = 1000000 * VECTOR(w)[i];

        if (v1 != v2) {
            // taking into account rounding up/down issues
            // e.g. for intermediate case no weights
            // 100000 == 99999
            if (v2+1 != v1) {
                return 0;
            }
        }
    }
    return 1;
}

int simple_test_case_no_weights_undirected() {

    igraph_t g;
    igraph_vector_t vector_edges, vector_expected_results, vector_actual_results;
    igraph_vector_t vector_expected_results_normalised;
    igraph_vector_init(&vector_actual_results, 3);

    igraph_real_t real_edges[] = {0,1 , 1,2};
    igraph_real_t real_expected_results[] = {0.333333, 0.5, 0.333333};
    igraph_real_t real_expected_results_normalised[] = {0.666666, 1, 0.666666};

    igraph_vector_view(&vector_edges, real_edges, sizeof(real_edges)/sizeof(igraph_real_t));
    igraph_create(&g, &vector_edges, /*number of vertices*/ 2, IGRAPH_DIRECTED);
    
    igraph_vector_view(&vector_expected_results, real_expected_results, 
                        sizeof(real_expected_results)/sizeof(igraph_real_t));
    
    igraph_vector_view(&vector_expected_results_normalised, 
                       real_expected_results_normalised, 
                       sizeof(real_expected_results_normalised)/sizeof(igraph_real_t));

    // NOT NORMALISED

    igraph_closeness_estimate(&g, &vector_actual_results /*store results here*/, 
                             /*calculating for all vectors in the graph*/ igraph_vss_all(), 
                             IGRAPH_ALL  /*graph is "undirected"*/, /*cutoff*/ 0 /*calculating exact centrality*/, 
                             NULL /*unweighted*/, /*not normalised*/ 0);

    if (igraph_check_vectors_equal(vector_expected_results, 
                                   vector_actual_results, 3) == 0) {
        igraph_vector_destroy(&vector_actual_results);
        igraph_destroy(&g);
        return 1;
    }

    // NORMALISED

    igraph_closeness_estimate(&g, &vector_actual_results /*store results here*/, 
                             /*calculating for all vectors in the graph*/ igraph_vss_all(), 
                             IGRAPH_ALL  /*graph is "undirected"*/, /*cutoff*/ 0 /*calculating exact centrality*/, 
                             NULL, /*normalised*/ 1);

    if (igraph_check_vectors_equal(vector_expected_results_normalised, 
                                   vector_actual_results, 3) == 0) {
        
        igraph_vector_destroy(&vector_actual_results);
        igraph_destroy(&g);
        return 1;
    }

    igraph_vector_destroy(&vector_actual_results);
    igraph_destroy(&g);

    return 0;
}

int simple_test_case_with_weights_undirected() {

    igraph_t g;
    igraph_vector_t vector_edges, vector_weights, vector_expected_results;
    igraph_vector_t vector_expected_results_normalised;
    igraph_vector_t vector_expected_results_cutoff;
    igraph_vector_t vector_actual_results;
    igraph_vector_init(&vector_actual_results, 3);

    igraph_real_t real_edges[] = {0,1 , 1,2};
    igraph_real_t real_weights[] = {3, 5};
    igraph_real_t real_expected_results[] = {0.090909, 0.125, 0.076923};
    igraph_real_t real_expected_results_normalised[] = {0.181818, 0.25, 0.153846};

    igraph_vector_view(&vector_edges, real_edges, sizeof(real_edges)/sizeof(igraph_real_t));
    igraph_create(&g, &vector_edges, /*number of vertices*/ 2, IGRAPH_DIRECTED);
    
    igraph_vector_view(&vector_expected_results, real_expected_results, 
                       sizeof(real_expected_results)/sizeof(igraph_real_t));
    
    igraph_vector_view(&vector_weights, real_weights, 
                       sizeof(real_weights)/sizeof(igraph_real_t));
    
    igraph_vector_view(&vector_expected_results_normalised, 
                       real_expected_results_normalised, 
                       sizeof(real_expected_results_normalised)/sizeof(igraph_real_t));
 
    // NOT NORMALISED

    igraph_closeness_estimate(&g, &vector_actual_results /*store results here*/, 
                             /*calculating for all vectors in the graph*/ igraph_vss_all(), 
                             IGRAPH_ALL  /*graph is "undirected"*/, /*cutoff*/ 0 /*calculating exact centrality*/, 
                             &vector_weights, /*not normalised*/ 0);

    if (igraph_check_vectors_equal(vector_expected_results, 
                                   vector_actual_results, 3) == 0) {
        
        igraph_vector_destroy(&vector_actual_results);
        igraph_destroy(&g);
        
        return 1;
    }

    // NORMALISED

    igraph_closeness_estimate(&g, &vector_actual_results /*store results here*/, 
                             /*calculating for all vectors in the graph*/ igraph_vss_all(), 
                             IGRAPH_ALL  /*graph is "undirected"*/, /*cutoff*/ 0 /*calculating exact centrality*/, 
                             &vector_weights, /*normalised*/ 1);

    if (igraph_check_vectors_equal(vector_expected_results_normalised, 
                                   vector_actual_results, 3) == 0) {
        igraph_vector_destroy(&vector_actual_results);
        igraph_destroy(&g);
        
        return 1;
    }

    igraph_vector_destroy(&vector_actual_results);
    igraph_destroy(&g);

    return 0;
}

int advanced_test_case_no_weights_undirected() {

    igraph_t g;
    igraph_vector_t vector_edges, vector_expected_results, vector_actual_results;
    igraph_vector_t vector_expected_results_normalised;
    igraph_vector_init(&vector_actual_results, 8);

    igraph_real_t real_edges[] = {1,0 , 0,5 , 5,6 , 5,4, 4,1 , 1,2 , 2,4 , 4,6 ,
                                 2,3 , 3,7 , 7,6 , 2,6};

    // note, denominatory calculated as (shortest dist)*(n-1) for no normalisation
    // normalisation excludes n-1 as part of the denominator
    igraph_real_t real_expected_results[] = {0.071428, 0.083333, 0.10, 0.071428,
                                             0.10, 0.083333, 0.10, 0.071428};
    igraph_real_t real_expected_results_normalised[] = {0.5, 0.583333, 0.70, 0.5,
                                                        0.70, 0.583333, 0.70, 0.5};    

    igraph_vector_view(&vector_edges, real_edges, sizeof(real_edges)/sizeof(igraph_real_t));
    igraph_create(&g, &vector_edges, /*number of vertices*/ 2, IGRAPH_DIRECTED);
    
    igraph_vector_view(&vector_expected_results, real_expected_results, 
                        sizeof(real_expected_results)/sizeof(igraph_real_t));
    
    igraph_vector_view(&vector_expected_results_normalised, 
                        real_expected_results_normalised, 
                        sizeof(real_expected_results_normalised)/sizeof(igraph_real_t));

    // NOT NORMALISED

    igraph_closeness_estimate(&g, &vector_actual_results /*store results here*/, 
                             /*calculating for all vectors in the graph*/ igraph_vss_all(), 
                             IGRAPH_ALL  /*graph is "undirected"*/, /*cutoff*/ 0 /*calculating exact centrality*/, 
                             NULL, /*not normalised*/ 0);

    if (igraph_check_vectors_equal(vector_expected_results, 
                                   vector_actual_results, 8) == 0) {
        
        igraph_vector_destroy(&vector_actual_results);
        igraph_destroy(&g);
        return 1;
    }

    // NORMALISED

    igraph_closeness_estimate(&g, &vector_actual_results /*store results here*/, 
                             /*calculating for all vectors in the graph*/ igraph_vss_all(), 
                             IGRAPH_ALL  /*graph is "undirected"*/, /*cutoff*/ 0 /*calculating exact centrality*/, 
                             NULL, /*normalised*/ 1);

    if (igraph_check_vectors_equal(vector_expected_results_normalised, 
                                   vector_actual_results, 3) == 0) {
        
        igraph_vector_destroy(&vector_actual_results);
        igraph_destroy(&g);
        
        return 1;
    }

    igraph_vector_destroy(&vector_actual_results);
    igraph_destroy(&g);

    return 0;
}

int advanced_test_case_with_weights() {

    igraph_t g;
    igraph_vector_t vector_edges, vector_weights, vector_expected_results_undirected;
    igraph_vector_t vector_expected_results_directed_out, vector_actual_results;
    igraph_vector_t vector_expected_results_directed_in;
    igraph_vector_init(&vector_actual_results, 8);

    igraph_real_t real_edges[] = {1,0 , 0,5 , 5,6 , 5,4, 4,1 , 1,2 , 2,4 , 4,6 ,
                                 2,3 , 3,7 , 7,6 , 6,2};
    
    igraph_real_t real_weights[] = {4, 9, 2, 2, 2, 3, 1, 1, 8, 7, 5, 5};
    
    igraph_real_t real_expected_results_undirected[] = {0.016949, 0.028571, 0.032258, 
                                             0.014084, 0.037037, 0.027027, 
                                             0.033333, 0.019230};
    igraph_real_t real_expected_results_directed_out[] = {0.008695, 0.017241, 0.019230, 
                                             0.007633, 0.016129, 0.016666, 
                                             0.011764, 0.010000};
    igraph_real_t real_expected_results_directed_in[] = {0.012820, 0.015873, 0.015873, 
                                             0.009803, 0.018867, 0.007518, 
                                             0.026315, 0.007518};


    igraph_vector_view(&vector_edges, real_edges, sizeof(real_edges)/sizeof(igraph_real_t));
    igraph_create(&g, &vector_edges, /*number of vertices*/ 2, IGRAPH_DIRECTED);
    
    igraph_vector_view(&vector_expected_results_undirected, 
                        real_expected_results_undirected, 
                        sizeof(real_expected_results_undirected)/sizeof(igraph_real_t));
    
    igraph_vector_view(&vector_expected_results_directed_out, 
                        real_expected_results_directed_out, 
                        sizeof(real_expected_results_undirected)/sizeof(igraph_real_t));
    
    igraph_vector_view(&vector_expected_results_directed_in, 
                        real_expected_results_directed_in, 
                        sizeof(real_expected_results_undirected)/sizeof(igraph_real_t));

    igraph_vector_view(&vector_weights, real_weights, sizeof(real_weights)/sizeof(igraph_real_t));

    // TEST FOR UNDIRECTED GRAPH

    igraph_closeness_estimate(&g, &vector_actual_results /*store results here*/, 
                             /*calculating for all vectors in the graph*/ igraph_vss_all(), 
                             IGRAPH_ALL  /*graph is "undirected"*/, /*cutoff*/ 0 /*calculating exact centrality*/, 
                             &vector_weights, /*not normalised*/ 0);

    if (igraph_check_vectors_equal(vector_expected_results_undirected, 
                                   vector_actual_results, 8) == 0) {
        
        igraph_vector_destroy(&vector_actual_results);
        igraph_destroy(&g);
        
        return 1;
    }

    // TEST FOR DIRECTED GRAPH
    // OUT means the min distance from the curr node to the other node

    igraph_closeness_estimate(&g, &vector_actual_results /*store results here*/, 
                             /*calculating for all vectors in the graph*/ igraph_vss_all(), 
                             IGRAPH_OUT  /*graph is "out directed"*/, /*cutoff*/ 0 /*calculating exact centrality*/, 
                             &vector_weights, /*not normalised*/ 0);

    if (igraph_check_vectors_equal(vector_expected_results_directed_out, 
                                   vector_actual_results, 8) == 0) {
        
        igraph_vector_destroy(&vector_actual_results);
        igraph_destroy(&g);
        
        return 1;
    }

    // IN means the min distance from a node to the curr node

    igraph_closeness_estimate(&g, &vector_actual_results /*store results here*/, 
                             /*calculating for all vectors in the graph*/ igraph_vss_all(), 
                             IGRAPH_IN  /*graph is "in directed"*/, /*cutoff*/ 0 /*calculating exact centrality*/, 
                             &vector_weights, /*not normalised*/ 0);

    if (igraph_check_vectors_equal(vector_expected_results_directed_in, 
                                   vector_actual_results, 8) == 0) {
        
        igraph_vector_destroy(&vector_actual_results);
        igraph_destroy(&g);
        
        return 1;
    }

    igraph_vector_destroy(&vector_actual_results);
    igraph_destroy(&g);

    return 0;
}

int disconnected_graph_test() {

    igraph_t g;
    igraph_vector_t vector_edges, vector_weights, vector_expected_results_unweighted;
    igraph_vector_t vector_actual_results, vector_expected_results_weighted;
    igraph_vector_init(&vector_actual_results, 10);

    // involves 2 disconnected subgraphs
    igraph_real_t real_edges[] = {0,1 , 1,2 , 2,3 , 1,3, 4,5 , 5,6 , 6,7 , 7,8 ,
                                 8,9 , 5,9};
    igraph_real_t real_expected_results_unweighted[] = {0.015384, 0.015873, 0.015625, 
                                             0.015625, 0.019607, 0.021276, 
                                             0.020833, 0.020408, 0.020408, 
                                             0.020833};
    
    igraph_real_t real_weights[] = {1, 2, 3, 7, 1, 1, 2, 4, 3, 12};
    igraph_real_t real_expected_results_weighted[] = {0.014285, 0.014705, 0.014705, 
                                             0.013513, 0.015151, 0.016129, 
                                             0.016666, 0.016666, 0.014705, 
                                             0.0125};

    igraph_vector_view(&vector_edges, real_edges, sizeof(real_edges)/sizeof(igraph_real_t));
    igraph_create(&g, &vector_edges, /*number of vertices*/ 2, IGRAPH_DIRECTED);
    
    igraph_vector_view(&vector_expected_results_unweighted, 
                        real_expected_results_unweighted, 
                        sizeof(real_expected_results_unweighted)/sizeof(igraph_real_t));
    
    igraph_vector_view(&vector_weights, real_weights, 
                        sizeof(real_weights)/sizeof(igraph_real_t)); 
    
    igraph_vector_view(&vector_expected_results_weighted, 
                        real_expected_results_weighted, 
                        sizeof(real_expected_results_weighted)/sizeof(igraph_real_t));

    // unweighted
    igraph_closeness_estimate(&g, &vector_actual_results /*store results here*/, 
                             /*calculating for all vectors in the graph*/ igraph_vss_all(), 
                             IGRAPH_ALL  /*graph is "undirected"*/, /*cutoff*/ 0 /*calculating exact centrality*/, 
                             NULL, /*not normalised*/ 0);

    if (igraph_check_vectors_equal(vector_expected_results_unweighted, 
                                   vector_actual_results, 10) == 0) {
        
        igraph_vector_destroy(&vector_actual_results);
        igraph_destroy(&g);
        
        return 1;
    }

    // weighted
    igraph_closeness_estimate(&g, &vector_actual_results /*store results here*/, 
                             /*calculating for all vectors in the graph*/ igraph_vss_all(), 
                             IGRAPH_ALL  /*graph is "undirected"*/, /*cutoff*/ 0 /*calculating exact centrality*/, 
                             &vector_weights, /*not normalised*/ 0);

    if (igraph_check_vectors_equal(vector_expected_results_weighted, 
                                   vector_actual_results, 10) == 0) {
        
        igraph_vector_destroy(&vector_actual_results);
        igraph_destroy(&g);
        
        return 1;
    }

    igraph_vector_destroy(&vector_actual_results);
    igraph_destroy(&g);

    return 0;

}

int linear_graph_test() {

    // testing a directed graph with a linear path, to check if the number of
    // vertices is used for the given path

    igraph_t g;
    igraph_vector_t vector_edges, vector_weights, vector_expected_results;
    igraph_vector_t vector_actual_results;
    igraph_vector_init(&vector_actual_results, 5);

    igraph_real_t real_edges[] = {0,1 , 1,2 , 2,3 , 3,4}; // 0-1-2-3-4
    igraph_real_t real_expected_results[] = {0.1, 0.090909, 0.076923, 0.0625, 0.05};
    
    igraph_vector_view(&vector_edges, real_edges, sizeof(real_edges)/sizeof(igraph_real_t));
    igraph_create(&g, &vector_edges, /*number of vertices*/ 2, IGRAPH_DIRECTED);
    
    igraph_vector_view(&vector_expected_results, real_expected_results, 
                       sizeof(real_expected_results)/sizeof(igraph_real_t));

    igraph_closeness_estimate(&g, &vector_actual_results /*store results here*/, 
                             /*calculating for all vectors in the graph*/ igraph_vss_all(), 
                             IGRAPH_OUT  /*graph is "out directed"*/, /*cutoff*/ 0 /*calculating exact centrality*/, 
                             NULL, /*not normalised*/ 0);

    if (igraph_check_vectors_equal(vector_expected_results, 
                                   vector_actual_results, 5) == 0) {
        
        igraph_vector_destroy(&vector_actual_results);
        igraph_destroy(&g);
        
        return 1;
    }

    igraph_vector_destroy(&vector_actual_results);
    igraph_destroy(&g);

    return 0;
}

int main(void) {

    if (simple_test_case_no_weights_undirected() == 1) {
        return 1;
    }
    if (simple_test_case_with_weights_undirected() == 1) {
        return 1;
    }
    if (advanced_test_case_no_weights_undirected() == 1) {
        return 1;
    }
    if (advanced_test_case_with_weights() == 1) {
        return 1;
    }
    if (disconnected_graph_test() == 1) {
        return 1;
    }
    if (linear_graph_test() == 1) {
        return 1;
    }
    
    return 0;
}