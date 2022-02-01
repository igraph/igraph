#include <igraph.h>
#include <stdio.h>

void print_vector(igraph_vector_t *vec){
    igraph_integer_t length = igraph_vector_size(vec);
    for(int i = 0; i < length; i++){
        printf(" %g", igraph_vector_e(vec, i));
    }
    printf("\n");
}

int main(){
    igraph_t g;
    igraph_vector_t degree, out_degree, in_degree;
    igraph_real_t assortativity;
    igraph_integer_t edge_count;

    /* Create graph */
    igraph_small(&g, 10, IGRAPH_DIRECTED, 
                 0, 1, 0, 2, 0, 3, 0, 5,
                 3, 6,
                 4, 5,
                 6, 4, 6, 7,
                 7, 8, 7, 9,
                 8, 9,
                 9, 0,
                 -1);
    
    edge_count = igraph_ecount(&g);

    /* Undirected graph */
    igraph_vector_init(&degree, 0);
    igraph_degree(&g, &degree, igraph_vss_all(), IGRAPH_ALL, /* consider self-loops */ IGRAPH_LOOPS);
    igraph_vector_sort(&degree);
    printf("Undirected graph degree sequence: ");
    print_vector(&degree);

    /* Assortativities for 5 different undirected graphs with identical degree sequences */
    for(int i = 0 ; i < 5 ; i++){
        igraph_rewire(&g, 4 * edge_count, IGRAPH_REWIRING_SIMPLE_LOOPS);
        igraph_assortativity(&g, &degree, NULL, &assortativity, /* ignore edge directions */ IGRAPH_UNDIRECTED);
        printf("Assortativity of undirected graph = %g\n", assortativity);
    }

    igraph_vector_destroy(&degree);
    printf("\n");

    /* Directed graph */
    igraph_vector_init(&out_degree, 0);
    igraph_degree(&g, &out_degree, igraph_vss_all(), IGRAPH_OUT, /* consider self-loops */ IGRAPH_LOOPS);
    igraph_vector_sort(&out_degree);
    printf("Directed graph degree sequence: ");
    print_vector(&out_degree);

    igraph_vector_init(&in_degree, 0);
    igraph_degree(&g, &in_degree, igraph_vss_all(), IGRAPH_IN, /* consider self-loops */ IGRAPH_LOOPS);
    igraph_vector_sort(&in_degree);
    printf("Directed graph degree sequence: ");
    print_vector(&in_degree);

    /* Assortativities for 5 different directed graphs with identical degree sequences */
    for(int i = 0 ; i < 5 ; i++){
        igraph_rewire(&g, 4 * edge_count, IGRAPH_REWIRING_SIMPLE_LOOPS);
        igraph_assortativity(&g, &out_degree, &in_degree, &assortativity, /* consider edge directions */ IGRAPH_DIRECTED);
        printf("Assortativity of undirected graph = %g\n", assortativity);
    }

    igraph_vector_destroy(&out_degree);
    igraph_vector_destroy(&in_degree);
    igraph_destroy(&g);
}
