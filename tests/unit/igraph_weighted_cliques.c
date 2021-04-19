
#include <igraph.h>
#include <stdlib.h>

#include "test_utilities.inc"

int compare_vectors(const void *p1, const void *p2) {
    igraph_vector_t *v1, *v2;
    long s1, s2, i;

    v1 = *((igraph_vector_t **) p1);
    v2 = *((igraph_vector_t **) p2);
    s1 = igraph_vector_size(v1);
    s2 = igraph_vector_size(v2);
    if (s1 < s2) {
        return -1;
    }
    if (s1 > s2) {
        return 1;
    }
    for (i = 0; i < s1; ++i) {
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
void canonicalize_list(igraph_vector_ptr_t *list) {
    long i, len;
    len = igraph_vector_ptr_size(list);
    for (i = 0; i < len; ++i) {
        igraph_vector_sort((igraph_vector_t *) VECTOR(*list)[i]);
    }
    qsort(&(VECTOR(*list)[0]), len, sizeof(void *), &compare_vectors);
}

/* Prints a clique vector along with its weight */
void print_weighted_clique(const igraph_vector_t *clique, const igraph_vector_t *vertex_weights) {
    long int i, n = igraph_vector_size(clique);
    igraph_real_t clique_weight = 0.0;
    for (i = 0; i < n; i++) {
        int v = VECTOR(*clique)[i];
        clique_weight += igraph_vector_e(vertex_weights, v);
        printf(" %d", v);
    }
    printf(" w=%.1f\n", clique_weight);
}

int main() {
    igraph_t graph;

    const igraph_integer_t n = 10; /* number of vertices in test graph */

    /* edges of the test graph */
    igraph_vector_t edges;
    igraph_real_t edge_data[] = {0., 1., 0., 6., 0., 7., 0., 8., 0., 9., 1., 2., 1., 3., 1., 4., 1.,
                                 6., 1., 7., 1., 8., 1., 9., 2., 3., 2., 5., 2., 6., 2., 7., 2., 9.,
                                 3., 5., 3., 6., 3., 7., 3., 9., 4., 5., 4., 6., 4., 7., 4., 9., 5.,
                                 8., 6., 7., 6., 8., 7., 8., 8., 9.
                                };

    /* vertex weights in test graph,
       note that current implementation only supports integer weights */
    igraph_vector_t vertex_weights;
    igraph_real_t vertex_weight_data[] = {3., 2., 3., 5., 2., 3., 1., 3., 5., 5.};

    igraph_vector_ptr_t result; /* result clique list */
    igraph_integer_t count; /* number of cliques found */

    igraph_real_t weighted_clique_no;

    int i;


    /* create graph */
    igraph_vector_init_copy(&edges, edge_data, (sizeof edge_data) / sizeof(igraph_real_t));
    igraph_create(&graph, &edges, n, /* directed= */ 0);

    /* set up vertex weight vector */
    igraph_vector_init_copy(&vertex_weights, vertex_weight_data, (sizeof vertex_weight_data) / sizeof(igraph_real_t));

    /* initialize result vector_ptr */
    igraph_vector_ptr_init(&result, 0);


    /* all weighed cliques above weight 6 */
    igraph_weighted_cliques(&graph, &vertex_weights, &result, 6, 0, /* maximal= */ 0);

    count = igraph_vector_ptr_size(&result);
    printf("%ld weighted cliques found above weight 6\n", (long) count);

    canonicalize_list(&result);
    for (i = 0; i < count; i++) {
        igraph_vector_t* v = (igraph_vector_t*) igraph_vector_ptr_e(&result, i);
        print_weighted_clique(v, &vertex_weights);
        igraph_vector_destroy(v);
        igraph_free(v);
    }


    /* all weighed cliques beteen weights 5 and 10 */
    igraph_weighted_cliques(&graph, &vertex_weights, &result, 5, 10, /* maximal= */ 0);

    count = igraph_vector_ptr_size(&result);
    printf("%ld weighted cliques found between weights 5 and 10\n", (long) count);

    canonicalize_list(&result);
    for (i = 0; i < count; i++) {
        igraph_vector_t* v = (igraph_vector_t*) igraph_vector_ptr_e(&result, i);
        print_weighted_clique(v, &vertex_weights);
        igraph_vector_destroy(v);
        igraph_free(v);
    }


    /* maximal weighed cliques above weight 7 */
    igraph_weighted_cliques(&graph, &vertex_weights, &result, 7, 0, /* maximal= */ 1);

    count = igraph_vector_ptr_size(&result);
    printf("%ld maximal weighted cliques found above weight 7\n", (long) count);

    canonicalize_list(&result);
    for (i = 0; i < count; i++) {
        igraph_vector_t* v = (igraph_vector_t*) igraph_vector_ptr_e(&result, i);
        print_weighted_clique(v, &vertex_weights);
        igraph_vector_destroy(v);
        igraph_free(v);
    }


    /* maximal weighed cliques beteen weights 5 and 10 */
    igraph_weighted_cliques(&graph, &vertex_weights, &result, 5, 10, /* maximal= */ 1);

    count = igraph_vector_ptr_size(&result);
    printf("%ld maximal weighted cliques found between weights 5 and 10\n", (long) count);

    canonicalize_list(&result);
    for (i = 0; i < count; i++) {
        igraph_vector_t* v = (igraph_vector_t*) igraph_vector_ptr_e(&result, i);
        print_weighted_clique(v, &vertex_weights);
        igraph_vector_destroy(v);
        igraph_free(v);
    }


    /* largest weight cliques */
    igraph_largest_weighted_cliques(&graph, &vertex_weights, &result);

    count = igraph_vector_ptr_size(&result);
    printf("%ld largest weight cliques found\n", (long) count);

    canonicalize_list(&result);
    for (i = 0; i < count; i++) {
        igraph_vector_t* v = (igraph_vector_t*) igraph_vector_ptr_e(&result, i);
        print_weighted_clique(v, &vertex_weights);
        igraph_vector_destroy(v);
        igraph_free(v);
    }

    igraph_weighted_clique_number(&graph, &vertex_weights, &weighted_clique_no);
    printf("weighted clique number: %.1f\n", weighted_clique_no);


    /* free data structures */
    igraph_vector_ptr_destroy(&result);
    igraph_vector_destroy(&vertex_weights);
    igraph_destroy(&graph);
    igraph_vector_destroy(&edges);

    VERIFY_FINALLY_STACK();

    return 0;
}
