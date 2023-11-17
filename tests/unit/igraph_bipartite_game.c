
#include <igraph.h>

#include "test_utilities.h"

void check_partitions(const igraph_t *g, const igraph_vector_bool_t *types, igraph_neimode_t mode) {
    igraph_integer_t m = igraph_ecount(g);

    for (igraph_integer_t i=0; i < m; i++) {
        switch (mode) {
        case IGRAPH_OUT:
            IGRAPH_ASSERT(VECTOR(*types)[IGRAPH_FROM(g, i)] == false);
            IGRAPH_ASSERT(VECTOR(*types)[IGRAPH_TO(g, i)] == true);
            break;
        case IGRAPH_IN:
            IGRAPH_ASSERT(VECTOR(*types)[IGRAPH_FROM(g, i)] == true);
            IGRAPH_ASSERT(VECTOR(*types)[IGRAPH_TO(g, i)] == false);
            break;
        case IGRAPH_ALL:
            IGRAPH_ASSERT(VECTOR(*types)[IGRAPH_FROM(g, i)] != VECTOR(*types)[IGRAPH_TO(g, i)]);
            break;
        }
    }
}

int main(void) {
    igraph_t graph;
    igraph_vector_bool_t types;
    igraph_bool_t bipartite;
    igraph_integer_t n1, n2, m;
    igraph_neimode_t modes[] = { IGRAPH_OUT, IGRAPH_IN, IGRAPH_ALL };

    igraph_rng_seed(igraph_rng_default(), 947);

    igraph_vector_bool_init(&types, 0);

    /* G(n,m) */

    /* undirected */

    n1 = 10; n2 = 20; m = 80;
    igraph_bipartite_game_gnm(&graph, &types,
                          n1, n2, m,
                          IGRAPH_UNDIRECTED, IGRAPH_ALL);

    igraph_is_bipartite(&graph, &bipartite, NULL);

    IGRAPH_ASSERT(bipartite);
    IGRAPH_ASSERT(! igraph_is_directed(&graph));
    IGRAPH_ASSERT(igraph_vcount(&graph) == n1 + n2);
    IGRAPH_ASSERT(igraph_ecount(&graph) == m);

    check_partitions(&graph, &types, IGRAPH_ALL);

    igraph_destroy(&graph);

    /* complete graph */

    n1 = 5; n2 = 6; m = 30;
    igraph_bipartite_game_gnm(&graph, &types,
                              n1, n2, m,
                              IGRAPH_UNDIRECTED, IGRAPH_ALL);

    igraph_is_bipartite(&graph, &bipartite, NULL);

    IGRAPH_ASSERT(bipartite);
    IGRAPH_ASSERT(! igraph_is_directed(&graph));
    IGRAPH_ASSERT(igraph_vcount(&graph) == n1 + n2);
    IGRAPH_ASSERT(igraph_ecount(&graph) == m);

    check_partitions(&graph, &types, IGRAPH_ALL);

    igraph_destroy(&graph);

    /* empty graph */

    n1 = 5; n2 = 6; m = 0;
    igraph_bipartite_game_gnm(&graph, &types,
                              n1, n2, m,
                              IGRAPH_UNDIRECTED, IGRAPH_ALL);

    igraph_is_bipartite(&graph, &bipartite, NULL);

    IGRAPH_ASSERT(bipartite);
    IGRAPH_ASSERT(! igraph_is_directed(&graph));
    IGRAPH_ASSERT(igraph_vcount(&graph) == n1 + n2);
    IGRAPH_ASSERT(igraph_ecount(&graph) == m);

    check_partitions(&graph, &types, IGRAPH_ALL);

    igraph_destroy(&graph);

    /* directed */

    for (int i=0; i < sizeof(modes) / sizeof(modes[0]); i++) {
        igraph_bipartite_game_gnm(&graph, &types,
                                  n1, n2, m,
                                  IGRAPH_DIRECTED, modes[i]);

        igraph_is_bipartite(&graph, &bipartite, NULL);

        IGRAPH_ASSERT(bipartite);
        IGRAPH_ASSERT(igraph_is_directed(&graph));
        IGRAPH_ASSERT(igraph_vcount(&graph) == n1 + n2);
        IGRAPH_ASSERT(igraph_ecount(&graph) == m);

        check_partitions(&graph, &types, modes[i]);

        igraph_destroy(&graph);
    }

    /* G(n,p) */

    /* undirected */

    n1 = 8; n2 = 15;
    igraph_bipartite_game_gnp(&graph, &types,
                          n1, n2, 0.8,
                          IGRAPH_UNDIRECTED, IGRAPH_ALL);

    igraph_is_bipartite(&graph, &bipartite, NULL);

    IGRAPH_ASSERT(bipartite);
    IGRAPH_ASSERT(! igraph_is_directed(&graph));
    IGRAPH_ASSERT(igraph_vcount(&graph) == n1 + n2);
    IGRAPH_ASSERT(igraph_ecount(&graph) > 0); /* 0 is exceedingly unlikely */

    igraph_destroy(&graph);

    /* directed */

    for (int i=0; i < sizeof(modes) / sizeof(modes[0]); i++) {
        igraph_bipartite_game_gnp(&graph, &types,
                                  n1, n2, 0.8,
                                  IGRAPH_DIRECTED, modes[i]);

        igraph_is_bipartite(&graph, &bipartite, NULL);

        IGRAPH_ASSERT(bipartite);
        IGRAPH_ASSERT(igraph_is_directed(&graph));
        IGRAPH_ASSERT(igraph_vcount(&graph) == n1 + n2);
        IGRAPH_ASSERT(igraph_ecount(&graph) > 0); /* 0 is exceedingly unlikely */

        check_partitions(&graph, &types, modes[i]);

        igraph_destroy(&graph);
    }

    igraph_vector_bool_destroy(&types);

    VERIFY_FINALLY_STACK();

    CHECK_ERROR(igraph_bipartite_game_gnm(&graph, NULL, 0, 10, 20, IGRAPH_DIRECTED, IGRAPH_ALL), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_bipartite_game_gnm(&graph, NULL, 10, 10, 201, IGRAPH_DIRECTED, IGRAPH_ALL), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_bipartite_game_gnm(&graph, NULL, -1, 10, 20, IGRAPH_DIRECTED, IGRAPH_ALL), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_bipartite_game_gnm(&graph, NULL, 10, -1, 20, IGRAPH_DIRECTED, IGRAPH_ALL), IGRAPH_EINVAL);

    CHECK_ERROR(igraph_bipartite_game_gnp(&graph, NULL, -1, 10, 0.1, IGRAPH_UNDIRECTED, IGRAPH_ALL), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_bipartite_game_gnp(&graph, NULL, 10, -1, 0.9, IGRAPH_UNDIRECTED, IGRAPH_ALL), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_bipartite_game_gnp(&graph, NULL, 10, 10, 1.1, IGRAPH_UNDIRECTED, IGRAPH_ALL), IGRAPH_EINVAL);

    VERIFY_FINALLY_STACK();

    return 0;
}
