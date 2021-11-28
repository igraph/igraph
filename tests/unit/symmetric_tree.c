
#include <igraph.h>
#include "test_utilities.inc"


int main() {

    // TODO: Add a test
    igraph_t g;
    igraph_vector_int_t v;
    igraph_vector_int_init(&v, 3);
    VECTOR(v)[0] = 3;
    VECTOR(v)[1] = 4;
    VECTOR(v)[2] = 5;


    IGRAPH_ASSERT(igraph_symmetric_tree(&g, &v, IGRAPH_TREE_OUT) == IGRAPH_SUCCESS);

    igraph_vector_int_init(&v, 4);
    VECTOR(v)[0] = 3;
    VECTOR(v)[1] = 4;
    VECTOR(v)[2] = 5;
    VECTOR(v)[3] = 6;

    // todo: add more tests

    return 0;
}