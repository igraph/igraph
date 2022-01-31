#include <igraph.h>
#include <stdio.h>

int main(){
    igraph_t g;
    igraph_vector_t types;
    igraph_real_t assortativity;

    /* Create directed graph */
    igraph_small(&g, 10, IGRAPH_DIRECTED, 
                 0, 1, 0, 2, 0, 3, 0, 5,
                 3, 6,
                 4, 5,
                 6, 4, 6, 7,
                 7, 8, 7, 9,
                 8, 9,
                 9, 0,
                 -1);

    igraph_vector_init(&types, 0);
    igraph_degree(&g, &types, igraph_vss_all(), IGRAPH_ALL, /* consider self-loops */ IGRAPH_LOOPS);

    igraph_assortativity(&g, &types, NULL, &assortativity, /* ignore edge directions */ IGRAPH_UNDIRECTED);
    printf("Assortativity of undirected graph = %g\n", assortativity);
    igraph_assortativity(&g, &types, NULL, &assortativity, /* consider edge directions */ IGRAPH_DIRECTED);
    printf("Assortativity of directed graph = %g\n", assortativity);
    
    igraph_vector_destroy(&types);
    igraph_destroy(&g);
}
