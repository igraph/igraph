/* -*- mode: C -*-  */

#include <igraph.h>

void print_vector(igraph_vector_t *v) {
    long int i, n = igraph_vector_size(v);
    for (i = 0; i < n; i++) {
        printf(" %li", (long int) VECTOR(*v)[i]);
    }
    printf("\n");
}

int main() {

    igraph_t g;
    igraph_real_t result;
    igraph_integer_t from, to;
    igraph_vector_t path, path_edge;

    igraph_barabasi_game(&g, 30, /*power=*/ 1, 30, 0, 0, /*A=*/ 1,
                         IGRAPH_DIRECTED, IGRAPH_BARABASI_BAG,
                         /*start_from=*/ 0);
    igraph_diameter(&g, &result, 0, 0, 0, 0, IGRAPH_UNDIRECTED, 1);

    /*   printf("Diameter: %li\n", (long int) result); */

    igraph_destroy(&g);

    igraph_ring(&g, 10, IGRAPH_DIRECTED, 0, 0);
    igraph_vector_init(&path, 0);
    igraph_vector_init(&path_edge, 0);
    igraph_diameter(&g, &result, &from, &to, &path, &path_edge, IGRAPH_DIRECTED, 1);
    printf("diameter: %li, from %li to %li\n", (long int) result,
           (long int) from, (long int) to);
    print_vector(&path);
    print_vector(&path_edge);

    igraph_vector_destroy(&path);
    igraph_vector_destroy(&path_edge);
    igraph_destroy(&g);

    return 0;
}