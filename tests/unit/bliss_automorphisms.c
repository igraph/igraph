
#include <igraph.h>

#include "test_utilities.inc"

#define TEST_GRAPH(name) \
    igraph_automorphisms(&graph, NULL, IGRAPH_BLISS_F, &info); \
    printf("%s: %s\n", name, info.group_size); \
    igraph_free(info.group_size); \
    igraph_destroy(&graph);

#define TEST_FAMOUS(name) \
    igraph_famous(&graph, name); \
    TEST_GRAPH(name);

int main() {
    igraph_t graph;
    igraph_bliss_info_t info;

    TEST_FAMOUS("Frucht");
    TEST_FAMOUS("Coxeter");
    TEST_FAMOUS("Petersen");
    TEST_FAMOUS("Meredith");

    igraph_full(&graph, 23, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    TEST_GRAPH("Complete 23");

    igraph_star(&graph, 17, IGRAPH_STAR_OUT, 0);
    TEST_GRAPH("Directed star 17");

    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    TEST_GRAPH("Null graph");

    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    TEST_GRAPH("Singleton graph");

    VERIFY_FINALLY_STACK();

    return 0;
}
