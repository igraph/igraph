
#include <igraph.h>

int main(void) {

    igraph_t graph;
    igraph_bool_t res;
    igraph_vector_int_t v;
    igraph_vector_int_init_int(&v, 3, 3, 4, 5);

    /* Initialize the library. */
    igraph_setup();

    /* Create a directed symmetric tree with 2 levels -
       3 children in first and 4 children in second level,
       5 children in third level
       with edges pointing towards the root. */
    igraph_symmetric_tree(&graph, &v, IGRAPH_TREE_IN);

    igraph_is_tree(&graph, &res, NULL, IGRAPH_IN);
    printf("Is it an in-tree? %s\n", res ? "Yes" : "No");

    igraph_is_tree(&graph, &res, NULL, IGRAPH_OUT);
    printf("Is it an out-tree? %s\n", res ? "Yes" : "No");

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&v);

    return 0;
}
