/* -*- mode: C -*-  */

#include <igraph.h>

void print_vector_int(igraph_vector_int_t *v) {
    igraph_integer_t i, n = igraph_vector_int_size(v);
    for (i = 0; i < n; i++) {
        printf(" %" IGRAPH_PRId, VECTOR(*v)[i]);
    }
    printf("\n");
}

int main(void) {

    igraph_t g;
    igraph_real_t result;
    igraph_integer_t from, to;
    igraph_vector_int_t path, path_edge;

    igraph_barabasi_game(&g, 30, /*power=*/ 1, 30, 0, 0, /*A=*/ 1,
                         IGRAPH_DIRECTED, IGRAPH_BARABASI_BAG,
                         /*start_from=*/ 0);
    igraph_diameter(&g, &result, 0, 0, 0, 0, IGRAPH_UNDIRECTED, 1);

    /*   printf("Diameter: %" IGRAPH_PRId "\n", (igraph_integer_t) result); */

    igraph_destroy(&g);

    igraph_ring(&g, 10, IGRAPH_DIRECTED, 0, 0);
    igraph_vector_int_init(&path, 0);
    igraph_vector_int_init(&path_edge, 0);
    igraph_diameter(&g, &result, &from, &to, &path, &path_edge, IGRAPH_DIRECTED, 1);
    printf(
        "diameter: %" IGRAPH_PRId ", from %" IGRAPH_PRId " to %" IGRAPH_PRId "\n",
        (igraph_integer_t) result, from, to
    );
    print_vector_int(&path);
    print_vector_int(&path_edge);

    igraph_vector_int_destroy(&path);
    igraph_vector_int_destroy(&path_edge);
    igraph_destroy(&g);

    return 0;
}
