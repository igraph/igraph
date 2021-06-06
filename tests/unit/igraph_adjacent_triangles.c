
#include <igraph.h>

#include "test_utilities.inc"

int main() {
    igraph_t g;
    igraph_vector_t res;
    igraph_vs_t vertices;

    igraph_vector_init(&res, 0);

    igraph_small(&g, 20, IGRAPH_DIRECTED,
                 15, 12, 12, 10, 15, 0, 11, 10, 2, 8, 8, 6, 13, 17, 10, 10, 17, 2, 14,
                 0, 16, 13, 14, 14, 0, 5, 6, 4, 0, 9, 0, 6, 10, 9, 16, 4, 14, 5, 17,
                 15, 14, 9, 17, 17, 1, 4, 10, 16, 7, 0, 11, 12, 6, 13, 2, 17, 4, 0, 0,
                 14, 4, 0, 6, 16, 16, 14, 13, 13, 12, 11, 3, 11, 11, 3, 6, 7, 4, 14,
                 10, 8, 13, 7, 14, 2, 5, 2, 0, 14, 3, 15, 5, 5, 7, 2, 14, 15, 5, 10,
                 10, 16, 7, 9, 14, 0, 15, 7, 13, 1, 15, 1, 4, 5, 4, 6, 16, 13, 6, 17,
                 8, 6, 9, 3, 8, 6, 6, 14, 11, 14, 6, 10, 10, 5, 1, 0, 16, 17, 9, 1, 5,
                 0, 5, 15, 8, 0, 0, 8, 5, 3, 9, 4, 13, 12, 11, 0, 11, 0, 10, 6, 4, 13,
                 8, 9, 11, 11, 3, 16, 1, 2, 16, 0, 9, 8, 3, 8, 8, 7, 12, 10, 9, 3, 13,
                 5, 3, 9, 6, 2, 11, 10, 1, 16, 0, 2, 10, 17, 16, 8, 11, 5, 13, 0, 19, 19,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 -1);

    igraph_vs_seq(&vertices, 0, igraph_vcount(&g) - 1);

    printf("\nDirected multi:\n");
    igraph_adjacent_triangles(&g, &res, igraph_vss_all());
    print_vector(&res);
    igraph_adjacent_triangles(&g, &res, vertices);
    print_vector(&res);

    printf("\nUndirected multi:\n");
    igraph_to_undirected(&g, IGRAPH_TO_UNDIRECTED_COLLAPSE, NULL);
    igraph_adjacent_triangles(&g, &res, igraph_vss_all());
    print_vector(&res);
    igraph_adjacent_triangles(&g, &res, vertices);
    print_vector(&res);

    printf("\nSimple:\n");
    igraph_simplify(&g, 1, 1, NULL);
    igraph_adjacent_triangles(&g, &res, igraph_vss_all());
    print_vector(&res);
    igraph_adjacent_triangles(&g, &res, vertices);
    print_vector(&res);

    igraph_vs_destroy(&vertices);

    igraph_destroy(&g);

    igraph_vector_destroy(&res);

    VERIFY_FINALLY_STACK();

    return 0;
}
