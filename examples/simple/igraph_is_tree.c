
#include <igraph.h>
#include <assert.h>

int main() {
    igraph_t g;
    igraph_bool_t res;
    igraph_integer_t root;

    /* the null graph is not a tree */
    igraph_empty(&g, 0, 0);

    igraph_is_tree(&g, &res, &root, IGRAPH_ALL);
    assert(! res);

    igraph_destroy(&g);

    /* the single-vertex graph is a tree */
    igraph_empty(&g, 1, 0);

    root = -1;
    igraph_is_tree(&g, &res, &root, IGRAPH_ALL);
    assert(res);
    assert(root == 0);

    igraph_destroy(&g);

    /*  4-cycle, not a tree */
    igraph_small(&g, 4, 0,
                 0,1, 1,2, 2,3, 3,0, -1);

    igraph_is_tree(&g, &res, &root, IGRAPH_ALL);
    assert(! res);

    igraph_destroy(&g);

    /* disconnected graph with the same number of edges a tree would have */
    igraph_small(&g, 4, 0,
                 0,1, 1,2, 0,2, 3,4, -1);

    igraph_is_tree(&g, &res, &root, IGRAPH_ALL);
    assert(! res);

    igraph_destroy(&g);

    /* 3-star, tree */
    igraph_small(&g, 4, 0,
                 0,1, 0,2, 0,3, -1);

    root = -1;
    igraph_is_tree(&g, &res, &root, IGRAPH_ALL);
    assert(res);
    assert(root == 0);

    igraph_destroy(&g);

    /* out-tree */
    igraph_small(&g, 4, 1,
                 0,1, 1,2, 1,3, -1);

    root = -1;
    igraph_is_tree(&g, &res, &root, IGRAPH_OUT);
    assert(res);
    assert(root == 0);

    igraph_is_tree(&g, &res, &root, IGRAPH_IN);
    assert(! res);

    root = -1;
    igraph_is_tree(&g, &res, &root, IGRAPH_ALL);
    assert(res);
    assert(root == 0);

    igraph_destroy(&g);

    /* in-tree */
    igraph_small(&g, 4, 1,
                 0,1, 2,1, 1,3, -1);

    root = -1;
    igraph_is_tree(&g, &res, &root, IGRAPH_IN);
    assert(res);
    assert(root == 3);

    igraph_is_tree(&g, &res, &root, IGRAPH_OUT);
    assert(! res);

    root = -1;
    igraph_is_tree(&g, &res, &root, IGRAPH_ALL);
    assert(res);
    assert(root == 0);

    igraph_destroy(&g);

    /* neither an in-tree, nor an out-ree, but still a tree when ignoring edge-directions */
    igraph_small(&g, 6, 1,
                 0,1, 1,2, 2,3, 4,1, 2,5, -1);

    root = -1;
    igraph_is_tree(&g, &res, &root, IGRAPH_ALL);
    assert(res);
    assert(root == 0);

    igraph_is_tree(&g, &res, &root, IGRAPH_IN);
    assert(! res);

    igraph_is_tree(&g, &res, &root, IGRAPH_OUT);
    assert(! res);

    igraph_destroy(&g);

    return 0;
}
