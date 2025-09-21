
#include <igraph.h>
#include <stdlib.h>

#include "test_utilities.h"

int compare_vectors(const igraph_vector_int_t *v1, const igraph_vector_int_t *v2) {
    igraph_int_t s1 = igraph_vector_int_size(v1);
    igraph_int_t s2 = igraph_vector_int_size(v2);
    if (s1 < s2) {
        return -1;
    }
    if (s1 > s2) {
        return 1;
    }
    for (igraph_int_t i = 0; i < s1; ++i) {
        if (VECTOR(*v1)[i] < VECTOR(*v2)[i]) {
            return -1;
        }
        if (VECTOR(*v1)[i] > VECTOR(*v2)[i]) {
            return 1;
        }
    }
    return 0;
}

/* Takes a pointer vector of vectors. Sorts each vector, then sorts the pointer vector */
void canonicalize_list(igraph_vector_int_list_t *list) {
    igraph_int_t len = igraph_vector_int_list_size(list);
    for (igraph_int_t i = 0; i < len; ++i) {
        igraph_vector_int_sort(igraph_vector_int_list_get_ptr(list, i));
    }
    igraph_vector_int_list_sort(list, &compare_vectors);
}

/* Prints a clique vector along with its weight */
void print_weighted_clique(const igraph_vector_int_t *clique, const igraph_vector_t *vertex_weights) {
    igraph_int_t i, n = igraph_vector_int_size(clique);
    igraph_real_t clique_weight = 0.0;
    for (i = 0; i < n; i++) {
        int v = VECTOR(*clique)[i];
        clique_weight += vertex_weights ? igraph_vector_get(vertex_weights, v) : 1;
        printf(" %d", v);
    }
    printf(" w=%g\n", clique_weight);
}

/* Prints a clique list and clears it */
void print_and_clear_weighted_clique_list(igraph_vector_int_list_t *cliques, const igraph_vector_t *vertex_weights) {
    igraph_int_t i, count;

    canonicalize_list(cliques);

    count = igraph_vector_int_list_size(cliques);
    for (i = 0; i < count; i++) {
        igraph_vector_int_t* v = igraph_vector_int_list_get_ptr(cliques, i);
        print_weighted_clique(v, vertex_weights);
    }

    igraph_vector_int_list_clear(cliques);
}

int main(void) {
    igraph_t graph;

    const igraph_int_t n = 10; /* number of vertices in test graph */

    /* edges of the test graph */
    igraph_vector_int_t edges;
    igraph_int_t edge_data[] = {0, 1, 0, 6, 0, 7, 0, 8, 0, 9, 1, 2, 1, 3, 1, 4, 1,
                                 6, 1, 7, 1, 8, 1, 9, 2, 3, 2, 5, 2, 6, 2, 7, 2, 9,
                                 3, 5, 3, 6, 3, 7, 3, 9, 4, 5, 4, 6, 4, 7, 4, 9, 5,
                                 8, 6, 7, 6, 8, 7, 8, 8, 9
                                };

    /* vertex weights in test graph,
       note that current implementation only supports integer weights */
    igraph_vector_t vertex_weights;
    igraph_real_t vertex_weight_data[] = {3., 2., 3., 5., 2., 3., 1., 3., 5., 5.};

    igraph_vector_int_list_t result; /* result clique list */
    igraph_int_t count; /* number of cliques found */

    igraph_real_t weighted_clique_no;


    /* create graph */
    igraph_vector_int_init_array(&edges, edge_data, (sizeof edge_data) / sizeof(edge_data[0]));
    igraph_create(&graph, &edges, n, IGRAPH_UNDIRECTED);

    /* set up vertex weight vector */
    igraph_vector_init_array(&vertex_weights, vertex_weight_data, (sizeof vertex_weight_data) / sizeof(vertex_weight_data[0]));

    /* initialize result vector_ptr */
    igraph_vector_int_list_init(&result, 0);


    /* all weighted cliques above weight 6 */
    igraph_weighted_cliques(&graph, &vertex_weights, &result, /* maximal= */ false, 6, 0, IGRAPH_UNLIMITED);

    count = igraph_vector_int_list_size(&result);
    printf("%" IGRAPH_PRId " weighted cliques found above weight 6\n", count);
    print_and_clear_weighted_clique_list(&result, &vertex_weights);


    /* all weighted cliques between weights 5 and 10 */
    igraph_weighted_cliques(&graph, &vertex_weights, &result, /* maximal= */ false, 5, 10, IGRAPH_UNLIMITED);

    count = igraph_vector_int_list_size(&result);
    printf("%" IGRAPH_PRId " weighted cliques found between weights 5 and 10\n", count);
    print_and_clear_weighted_clique_list(&result, &vertex_weights);


    /* maximal weighted cliques above weight 7 */
    igraph_weighted_cliques(&graph, &vertex_weights, &result, /* maximal= */ true, 7, 0, IGRAPH_UNLIMITED);

    count = igraph_vector_int_list_size(&result);
    printf("%" IGRAPH_PRId " maximal weighted cliques found above weight 7\n", count);
    print_and_clear_weighted_clique_list(&result, &vertex_weights);


    /* maximal weighed cliques beteen weights 5 and 10 */
    igraph_weighted_cliques(&graph, &vertex_weights, &result, /* maximal= */ true, 5, 10, IGRAPH_UNLIMITED);

    count = igraph_vector_int_list_size(&result);
    printf("%" IGRAPH_PRId " maximal weighted cliques found between weights 5 and 10\n", count);
    print_and_clear_weighted_clique_list(&result, &vertex_weights);


    /* largest weight cliques */
    igraph_largest_weighted_cliques(&graph, &vertex_weights, &result);

    count = igraph_vector_int_list_size(&result);
    printf("%" IGRAPH_PRId " largest weight cliques found\n", count);
    print_and_clear_weighted_clique_list(&result, &vertex_weights);

    igraph_weighted_clique_number(&graph, &vertex_weights, &weighted_clique_no);
    printf("weighted clique number: %g\n", weighted_clique_no);


    /* Test unweighted fallback: */
    printf("\nUnweighted case:\n");

    /* test fallback to unweighted variants: all cliques */
    igraph_weighted_cliques(&graph, NULL, &result, /* maximal= */ false, 4, 5, IGRAPH_UNLIMITED);

    count = igraph_vector_int_list_size(&result);
    printf("%" IGRAPH_PRId " unweighted cliques found between sizes 4 and 5\n", count);
    print_and_clear_weighted_clique_list(&result, 0);

    /* test fallback to unweighted variants: maximal cliques */
    igraph_weighted_cliques(&graph, NULL, &result, /* maximal= */ true, 4, 5, IGRAPH_UNLIMITED);

    count = igraph_vector_int_list_size(&result);
    printf("%" IGRAPH_PRId " unweighted maximal cliques found between sizes 4 and 5\n", count);
    print_and_clear_weighted_clique_list(&result, 0);

    /* test fallback to unweighted variants: largest cliques */
    igraph_largest_weighted_cliques(&graph, NULL, &result);

    count = igraph_vector_int_list_size(&result);
    printf("%" IGRAPH_PRId " largest unweighted cliques found\n", count);
    print_and_clear_weighted_clique_list(&result, 0);

    /* test fallback to unweighted variants: clique number */
    igraph_weighted_clique_number(&graph, NULL, &weighted_clique_no);
    printf("unweighted clique number: %g\n", weighted_clique_no);

    /* Here we test unit weights, which should give identical results to the unweighted case: */
    printf("\nUnit weights:\n");

    igraph_vector_fill(&vertex_weights, 1);

    /* test unit weights: all cliques */
    igraph_weighted_cliques(&graph, &vertex_weights, &result, /* maximal= */ false, 4, 5, IGRAPH_UNLIMITED);

    count = igraph_vector_int_list_size(&result);
    printf("%" IGRAPH_PRId " cliques with unit weights found between sizes 4 and 5\n", count);
    print_and_clear_weighted_clique_list(&result, 0);

    /* test unit weights: maximal cliques */
    igraph_weighted_cliques(&graph, &vertex_weights, &result, /* maximal= */ true, 4, 5, IGRAPH_UNLIMITED);

    count = igraph_vector_int_list_size(&result);
    printf("%" IGRAPH_PRId " maximal cliques with unit weights between sizes 4 and 5\n", count);
    print_and_clear_weighted_clique_list(&result, 0);

    /* test unit weights: largest cliques */
    igraph_largest_weighted_cliques(&graph, &vertex_weights, &result);

    count = igraph_vector_int_list_size(&result);
    printf("%" IGRAPH_PRId " largest cliques with unit weights found\n", count);
    print_and_clear_weighted_clique_list(&result, 0);

    /* test unit weights: clique number */
    igraph_weighted_clique_number(&graph, NULL, &weighted_clique_no);
    printf("clique number with unit weights: %g\n", weighted_clique_no);


    /* free data structures */
    igraph_vector_int_list_destroy(&result);
    igraph_vector_destroy(&vertex_weights);
    igraph_destroy(&graph);
    igraph_vector_int_destroy(&edges);

    /* additional small examples */

    printf("\nP_2 graph with weights (5, 5):\n");
    igraph_small(&graph, 2, IGRAPH_UNDIRECTED,
                 0,1,
                 -1);

    igraph_vector_init_int(&vertex_weights, 2,
                           5, 5);

    igraph_weighted_clique_number(&graph, &vertex_weights, &weighted_clique_no);
    printf("weighted clique number: %g\n", weighted_clique_no);

    printf("\nP_2 graph with weights (5, 4):\n");
    VECTOR(vertex_weights)[1] = 4;

    igraph_weighted_clique_number(&graph, &vertex_weights, &weighted_clique_no);
    printf("weighted clique number: %g\n", weighted_clique_no);

    igraph_vector_destroy(&vertex_weights);
    igraph_destroy(&graph);

    printf("\nK_3 graph with weights (3, 3, 3):\n");
    igraph_full(&graph, 3, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

    igraph_vector_init_int(&vertex_weights, 3,
                           3, 3, 3);

    igraph_weighted_clique_number(&graph, &vertex_weights, &weighted_clique_no);
    printf("weighted clique number: %g\n", weighted_clique_no);

    printf("\nK_3 graph with weights (3, 4, 3):\n");
    VECTOR(vertex_weights)[1] = 4;

    igraph_weighted_clique_number(&graph, &vertex_weights, &weighted_clique_no);
    printf("weighted clique number: %g\n", weighted_clique_no);

    igraph_vector_destroy(&vertex_weights);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
