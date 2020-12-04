
#include <igraph.h>
#include <stdlib.h>

#include "test_utilities.inc"

struct userdata {
    int i;
    igraph_vector_ptr_t *list;
};

int compare_vectors(const void *p1, const void *p2) {
    igraph_vector_t *v1, *v2;
    long s1, s2, i;

    v1 = *((igraph_vector_t **) p1);
    v2 = *((igraph_vector_t **) p2);
    s1 = igraph_vector_size(v1);
    s2 = igraph_vector_size(v2);
    if (s1 < s2) {
        return -1;
    }
    if (s1 > s2) {
        return 1;
    }
    for (i = 0; i < s1; ++i) {
        if (VECTOR(*v1)[i] < VECTOR(*v2)[i]) {
            return -1;
        }
        if (VECTOR(*v1)[i] > VECTOR(*v2)[i]) {
            return 1;
        }
    }
    return 0;
}


igraph_bool_t handler(igraph_vector_t *clique, void *arg) {
    struct userdata *ud;
    igraph_bool_t cont;

    ud = (struct userdata *) arg;
    cont = 1; /* true */

    if (compare_vectors(&clique, &(VECTOR(*(ud->list))[ud->i])) != 0) {
        printf("igraph_maximal_cliques() and igraph_maximal_cliques_callback() give different results.\n");
        cont = 0; /* false */
    }

    igraph_vector_destroy(clique);
    igraph_free(clique);

    ud->i += 1;

    return cont;
}


igraph_bool_t handler_stop(igraph_vector_t *clique, void *arg) {
    /* Stop search as soon as a 3-clique is found. */
    /* Since there are two 3-cliques in the test graph, this will stop the search before it is complete. */
    if (igraph_vector_size(clique) == 3) {
        igraph_vector_destroy(clique);
        igraph_free(clique);

        return 0;    /* false */
    }

    igraph_vector_destroy(clique);
    igraph_free(clique);

    return 1 /* true */;
}


int main() {
    igraph_t graph;
    igraph_vector_ptr_t list;
    struct userdata ud;

    igraph_small(&graph, 6, 0,
                 1, 2, 2, 3, 3, 4, 4, 5, 5, 2, 2, 4,
                 -1);

    igraph_vector_ptr_init(&list, 0);
    igraph_maximal_cliques(&graph, &list, 0, 0);

    ud.i = 0;
    ud.list = &list;

    /* Check that the callback function finds the same cliques as igraph_maximal_cliques() */
    igraph_maximal_cliques_callback(&graph, &handler, (void *) &ud, 0, 0);

    /* Check that the search can be stopped correctly */
    igraph_maximal_cliques_callback(&graph, &handler_stop, NULL, 0, 0);

    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&list, igraph_vector_destroy);
    igraph_vector_ptr_destroy_all(&list);

    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
