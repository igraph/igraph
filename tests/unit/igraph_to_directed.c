
#include <igraph.h>

#include "test_utilities.inc"

int main() {
    igraph_t ug;
    igraph_t dg;

    igraph_small(&ug, 6, /* directed= */ 0,
                 2,0, 1,4, 3,2, 2,3, 0,4, 5,0, 2,3, 5,3, 2,5, 0,1,
                 -1);

    igraph_copy(&dg, &ug);
    printf("\nARBITRARY:\n");
    igraph_to_directed(&dg, IGRAPH_TO_DIRECTED_ARBITRARY);
    print_graph(&dg);
    igraph_destroy(&dg);

    igraph_copy(&dg, &ug);
    printf("\nMUTUAL:\n");
    igraph_to_directed(&dg, IGRAPH_TO_DIRECTED_MUTUAL);
    print_graph(&dg);
    igraph_destroy(&dg);

    igraph_copy(&dg, &ug);
    printf("\nACYCLIC:\n");
    igraph_to_directed(&dg, IGRAPH_TO_DIRECTED_ACYCLIC);
    print_graph(&dg);
    igraph_destroy(&dg);

    igraph_copy(&dg, &ug);
    printf("\nRANDOM (edge count only):\n");
    igraph_to_directed(&dg, IGRAPH_TO_DIRECTED_RANDOM);
    printf("%" IGRAPH_PRId "\n", igraph_ecount(&dg));
    igraph_destroy(&dg);

    igraph_destroy(&ug);

    VERIFY_FINALLY_STACK();

    return 0;
}
