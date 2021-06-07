
#include <igraph.h>

#include "test_utilities.inc"

int main() {

    igraph_t g;
    igraph_real_t global, global2, global3;

    /* Small graphs */

    printf("Null graph: ");
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    igraph_transitivity_undirected(&g, &global, IGRAPH_TRANSITIVITY_NAN);
    print_real(stdout, global, "%g");
    printf("\n");
    igraph_destroy(&g);

    printf("\nSingleton graph: ");
    igraph_empty(&g, 1, IGRAPH_UNDIRECTED);
    igraph_transitivity_undirected(&g, &global, IGRAPH_TRANSITIVITY_NAN);
    print_real(stdout, global, "%g");
    printf("\n");
    igraph_destroy(&g);

    printf("\nTwo connected vertices: ");
    igraph_small(&g, 2, IGRAPH_UNDIRECTED, 0,1, -1);
    igraph_transitivity_undirected(&g, &global, IGRAPH_TRANSITIVITY_NAN);
    print_real(stdout, global, "%g");
    printf("\n");
    igraph_destroy(&g);

    printf("\nTriangle: ");
    igraph_full(&g, 3, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_transitivity_undirected(&g, &global, IGRAPH_TRANSITIVITY_NAN);
    print_real(stdout, global, "%g");
    printf("\n");
    igraph_destroy(&g);

    printf("\nTwo-star: ");
    igraph_small(&g, 3, IGRAPH_UNDIRECTED, 0,2, 0,1, -1);
    igraph_transitivity_undirected(&g, &global, IGRAPH_TRANSITIVITY_NAN);
    print_real(stdout, global, "%g");
    printf("\n");
    igraph_destroy(&g);

    printf("\nZachary karate club: ");
    igraph_famous(&g, "Zachary");
    igraph_transitivity_undirected(&g, &global, IGRAPH_TRANSITIVITY_NAN);
    print_real(stdout, global, "%g");
    printf("\n");
    igraph_destroy(&g);

    printf("\nDirected and multigraphs:\n");

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

    printf("\nDirected multi: ");
    igraph_transitivity_undirected(&g, &global, IGRAPH_TRANSITIVITY_NAN);
    print_real(stdout, global, "%.10g");
    printf("\n");

    printf("Undirected multi: ");
    igraph_to_undirected(&g, IGRAPH_TO_UNDIRECTED_COLLAPSE, NULL);
    igraph_transitivity_undirected(&g, &global2, IGRAPH_TRANSITIVITY_NAN);
    print_real(stdout, global2, "%.10g");
    printf("\n");

    printf("Simple: ");
    igraph_simplify(&g, 1, 1, NULL);
    igraph_transitivity_undirected(&g, &global3, IGRAPH_TRANSITIVITY_NAN);
    print_real(stdout, global3, "%.10g");
    printf("\n");

    IGRAPH_ASSERT(global == global2);
    IGRAPH_ASSERT(global == global3);

    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
