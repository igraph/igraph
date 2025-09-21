
#include <igraph.h>

int main(void) {
    igraph_t tree;
    igraph_vector_t eccentricity;
    igraph_bool_t is_tree;

    /* Initialize the library. */
    igraph_setup();

    /* Create a Bethe lattice with 5 levels, i.e. height 4. */
    igraph_regular_tree(&tree, 4, 3, IGRAPH_TREE_UNDIRECTED);

    /* Bethe lattices are trees. */
    igraph_is_tree(&tree, &is_tree, NULL, IGRAPH_ALL);
    printf("Is it a tree? %s\n", is_tree ? "Yes." : "No.");

    /* Compute and print eccentricities. The root is the most central. */
    igraph_vector_init(&eccentricity, 0);
    igraph_eccentricity(&tree, NULL, &eccentricity, igraph_vss_all(), IGRAPH_ALL);
    printf("Vertex eccentricities:\n");
    igraph_vector_print(&eccentricity);
    igraph_vector_destroy(&eccentricity);

    /* Clean up. */
    igraph_destroy(&tree);

    return 0;
}
