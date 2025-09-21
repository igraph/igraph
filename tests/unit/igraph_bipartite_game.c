/*
   igraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <igraph.h>

#include "test_utilities.h"

void check_partitions(
    const igraph_t *graph,
    const igraph_vector_bool_t *types,
    igraph_neimode_t mode
) {
    igraph_int_t m = igraph_ecount(graph);

    IGRAPH_ASSERT(igraph_vector_bool_size(types) == igraph_vcount(graph));

    for (igraph_int_t i=0; i < m; i++) {
        switch (mode) {
        case IGRAPH_OUT:
            IGRAPH_ASSERT(VECTOR(*types)[IGRAPH_FROM(graph, i)] == false);
            IGRAPH_ASSERT(VECTOR(*types)[IGRAPH_TO(graph, i)] == true);
            break;
        case IGRAPH_IN:
            IGRAPH_ASSERT(VECTOR(*types)[IGRAPH_FROM(graph, i)] == true);
            IGRAPH_ASSERT(VECTOR(*types)[IGRAPH_TO(graph, i)] == false);
            break;
        case IGRAPH_ALL:
            IGRAPH_ASSERT(VECTOR(*types)[IGRAPH_FROM(graph, i)] != VECTOR(*types)[IGRAPH_TO(graph, i)]);
            break;
        }
    }
}

void check_gnm(
    igraph_int_t n1, igraph_int_t n2, igraph_int_t m,
    igraph_bool_t directed, igraph_neimode_t mode, igraph_bool_t multi
) {
    igraph_t graph;
    igraph_vector_bool_t types;
    igraph_bool_t has_loop, has_multi;
    igraph_bool_t bipartite;

    if (! directed) {
        mode = IGRAPH_ALL;
    }

    igraph_vector_bool_init(&types, 0);
    igraph_bipartite_game_gnm(&graph, &types, n1, n2, m, directed, mode, multi ? IGRAPH_MULTI_SW : IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);

    /* check correct vertex and edge count, directedness */
    IGRAPH_ASSERT(igraph_is_directed(&graph) == directed);
    IGRAPH_ASSERT(igraph_vcount(&graph) == n1 + n2);
    IGRAPH_ASSERT(igraph_ecount(&graph) == m);

    /* bipartite graphs do not have self-loops */
    igraph_has_loop(&graph, &has_loop);
    IGRAPH_ASSERT(! has_loop);

    /* no multi-edges unless explicitly allowed */
    igraph_has_multiple(&graph, &has_multi);
    if (! multi) {
        IGRAPH_ASSERT(! has_multi);
    }

    /* redundant with the next check, but also tests is_bipartite() */
    igraph_is_bipartite(&graph, &bipartite, NULL);
    IGRAPH_ASSERT(bipartite);

    check_partitions(&graph, &types, mode);

    igraph_destroy(&graph);
    igraph_vector_bool_destroy(&types);
}

void check_iea(
    igraph_int_t n1, igraph_int_t n2, igraph_int_t m,
    igraph_bool_t directed, igraph_neimode_t mode
) {
    igraph_t graph;
    igraph_vector_bool_t types;
    igraph_bool_t has_loop;
    igraph_bool_t bipartite;

    if (! directed) {
        mode = IGRAPH_ALL;
    }

    igraph_vector_bool_init(&types, 0);
    igraph_bipartite_iea_game(&graph, &types, n1, n2, m, directed, mode);

    /* check correct vertex and edge count, directedness */
    IGRAPH_ASSERT(igraph_is_directed(&graph) == directed);
    IGRAPH_ASSERT(igraph_vcount(&graph) == n1 + n2);
    IGRAPH_ASSERT(igraph_ecount(&graph) == m);

    /* bipartite graphs do not have self-loops */
    igraph_has_loop(&graph, &has_loop);
    IGRAPH_ASSERT(! has_loop);

    /* redundant with the next check, but also tests is_bipartite() */
    igraph_is_bipartite(&graph, &bipartite, NULL);
    IGRAPH_ASSERT(bipartite);

    check_partitions(&graph, &types, mode);

    igraph_destroy(&graph);
    igraph_vector_bool_destroy(&types);
}

void check_gnp(
    igraph_int_t n1, igraph_int_t n2, igraph_real_t p,
    igraph_bool_t directed, igraph_neimode_t mode, igraph_bool_t multiple,
    igraph_bool_t edge_labeled,
    igraph_bool_t assume_edges
) {
    igraph_t graph;
    igraph_vector_bool_t types;
    igraph_bool_t has_loop, has_multi;
    igraph_bool_t bipartite;

    if (! directed) {
        mode = IGRAPH_ALL;
    }

    igraph_vector_bool_init(&types, 0);
    igraph_bipartite_game_gnp(&graph, &types, n1, n2, p, directed, mode, multiple ? IGRAPH_MULTI_SW : IGRAPH_SIMPLE_SW, edge_labeled);

    /* check correct vertex and edge count, directedness */
    IGRAPH_ASSERT(igraph_is_directed(&graph) == directed);
    IGRAPH_ASSERT(igraph_vcount(&graph) == n1 + n2);
    if (assume_edges) {
        /* with most parameter values, having no edges is exceedingly unlikely */
        IGRAPH_ASSERT(igraph_ecount(&graph) > 0);
    }
    if (p == 1) {
        /* complete graph */
        IGRAPH_ASSERT(igraph_ecount(&graph) == (directed && mode == IGRAPH_ALL) ? 2*n1*n2 : n1*n2);
    } else if (p == 0) {
        /* empty graph */
        IGRAPH_ASSERT(igraph_ecount(&graph) == 0);
    }

    /* bipartite graphs do not have self-loops */
    igraph_has_loop(&graph, &has_loop);
    IGRAPH_ASSERT(! has_loop);

    /* no multi-edges unless explicitly allowed */
    igraph_has_multiple(&graph, &has_multi);
    if (!multiple) IGRAPH_ASSERT(! has_multi);

    /* redundant with the next check, but also tests is_bipartite() */
    igraph_is_bipartite(&graph, &bipartite, NULL);
    IGRAPH_ASSERT(bipartite);

    check_partitions(&graph, &types, mode);

    igraph_destroy(&graph);
    igraph_vector_bool_destroy(&types);
}

int main(void) {
    igraph_t graph;
    igraph_neimode_t modes[] = { IGRAPH_OUT, IGRAPH_IN, IGRAPH_ALL };

    igraph_rng_seed(igraph_rng_default(), 947);

    /* G(n,m) */

    /* UNDIRECTED */

    /* null graph */
    check_gnm(0, 0, 0, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_NO_MULTIPLE);
    check_gnm(0, 0, 0, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_MULTIPLE);

    /* empty partition */
    check_gnm(0, 2, 0, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_NO_MULTIPLE);
    check_gnm(3, 0, 0, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_MULTIPLE);

    /* empty graph */
    check_gnm(5, 6, 0, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_NO_MULTIPLE);
    check_gnm(6, 5, 0, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_MULTIPLE);

    /* arbitrary graph */
    check_gnm(10, 20, 80, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_NO_MULTIPLE);
    check_gnm(10, 20, 80, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_MULTIPLE);
    check_gnm(20, 10, 80, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_MULTIPLE);
    check_gnm(8, 12, 150, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_MULTIPLE);
    check_gnm(12, 8, 150, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_MULTIPLE);

    /* complete graph */
    check_gnm(5, 6, 5*6, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_NO_MULTIPLE);
    check_gnm(6, 5, 5*6, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_NO_MULTIPLE);

    /* DIRECTED */

    for (size_t i=0; i < sizeof(modes) / sizeof(modes[0]); i++) {
        /* null graph */
        check_gnm(0, 0, 0, IGRAPH_DIRECTED, modes[i], IGRAPH_NO_MULTIPLE);
        check_gnm(0, 0, 0, IGRAPH_DIRECTED, modes[i], IGRAPH_MULTIPLE);

        /* empty partition */
        check_gnm(3, 0, 0, IGRAPH_DIRECTED, IGRAPH_ALL, IGRAPH_NO_MULTIPLE);
        check_gnm(0, 4, 0, IGRAPH_DIRECTED, IGRAPH_ALL, IGRAPH_MULTIPLE);

        /* empty graph */
        check_gnm(2, 4, 0, IGRAPH_DIRECTED, modes[i], IGRAPH_NO_MULTIPLE);
        check_gnm(5, 3, 0, IGRAPH_DIRECTED, modes[i], IGRAPH_MULTIPLE);

        /* arbitrary graph */

        check_gnm(4, 7, 15, IGRAPH_DIRECTED, modes[i], IGRAPH_NO_MULTIPLE);
        check_gnm(4, 7, 15, IGRAPH_DIRECTED, modes[i], IGRAPH_MULTIPLE);

        check_gnm(7, 4, 15, IGRAPH_DIRECTED, modes[i], IGRAPH_NO_MULTIPLE);
        check_gnm(7, 4, 15, IGRAPH_DIRECTED, modes[i], IGRAPH_MULTIPLE);

        /* many edges */
        check_gnm(4, 3, 25, IGRAPH_DIRECTED, modes[i], IGRAPH_MULTIPLE);
    }

    /* complete graph */
    check_gnm(3, 2, 6, IGRAPH_DIRECTED, IGRAPH_IN, IGRAPH_NO_MULTIPLE);
    check_gnm(3, 2, 12, IGRAPH_DIRECTED, IGRAPH_ALL, IGRAPH_NO_MULTIPLE);

    VERIFY_FINALLY_STACK();

    /* IEA */

    /* UNDIRECTED */

    /* null graph */
    check_iea(0, 0, 0, IGRAPH_UNDIRECTED, IGRAPH_ALL);

    /* empty partition */
    check_iea(2, 0, 0, IGRAPH_UNDIRECTED, IGRAPH_ALL);

    /* empty graph */
    check_iea(5, 6, 0, IGRAPH_UNDIRECTED, IGRAPH_ALL);

    /* arbitrary graph */
    check_iea(20, 10, 80, IGRAPH_UNDIRECTED, IGRAPH_ALL);
    check_iea(8, 12, 150, IGRAPH_UNDIRECTED, IGRAPH_ALL);

    /* DIRECTED */

    for (size_t i=0; i < sizeof(modes) / sizeof(modes[0]); i++) {
        /* null graph */
        check_iea(0, 0, 0, IGRAPH_DIRECTED, modes[i]);

        /* empty partition */
        check_iea(0, 2, 0, IGRAPH_DIRECTED, modes[i]);

        /* empty graph */
        check_iea(2, 4, 0, IGRAPH_DIRECTED, modes[i]);

        /* arbitrary graph */
        check_iea(4, 7, 15, IGRAPH_DIRECTED, modes[i]);

        /* many edges */
        check_iea(4, 3, 25, IGRAPH_DIRECTED, modes[i]);
    }

    VERIFY_FINALLY_STACK();

    /* G(n,p) */

    /* UNDIRECTED */

    check_gnp(0, 0, 0, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_NO_MULTIPLE, IGRAPH_EDGE_UNLABELED, false);
    check_gnp(2, 3, 0, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_NO_MULTIPLE, IGRAPH_EDGE_UNLABELED, false);
    check_gnp(0, 0, 0, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_MULTIPLE, IGRAPH_EDGE_UNLABELED, false);
    check_gnp(2, 3, 0, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_MULTIPLE, IGRAPH_EDGE_UNLABELED, false);
    check_gnp(0, 0, 0, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_MULTIPLE, IGRAPH_EDGE_LABELED, false);
    check_gnp(2, 3, 0, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_MULTIPLE, IGRAPH_EDGE_LABELED, false);

    check_gnp(8, 15, 0.8, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_NO_MULTIPLE, IGRAPH_EDGE_UNLABELED, true);
    check_gnp(6, 3, 1, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_NO_MULTIPLE, IGRAPH_EDGE_UNLABELED, true);
    check_gnp(8, 15, 0.8, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_MULTIPLE, IGRAPH_EDGE_UNLABELED, true);
    check_gnp(6, 3, 1, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_MULTIPLE, IGRAPH_EDGE_UNLABELED, true);
    check_gnp(8, 15, 0.8, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_MULTIPLE, IGRAPH_EDGE_LABELED, true);
    check_gnp(6, 3, 1, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_MULTIPLE, IGRAPH_EDGE_LABELED, true);

    /* DIRECTED */

    for (size_t i=0; i < sizeof(modes) / sizeof(modes[0]); i++) {
        check_gnp(0, 0, 0, IGRAPH_DIRECTED, modes[i], IGRAPH_NO_MULTIPLE, IGRAPH_EDGE_UNLABELED, false);
        check_gnp(2, 3, 0, IGRAPH_DIRECTED, modes[i], IGRAPH_NO_MULTIPLE, IGRAPH_EDGE_UNLABELED, false);
        check_gnp(0, 0, 0, IGRAPH_DIRECTED, modes[i], IGRAPH_MULTIPLE, IGRAPH_EDGE_UNLABELED, false);
        check_gnp(2, 3, 0, IGRAPH_DIRECTED, modes[i], IGRAPH_MULTIPLE, IGRAPH_EDGE_UNLABELED, false);
        check_gnp(0, 0, 0, IGRAPH_DIRECTED, modes[i], IGRAPH_MULTIPLE, IGRAPH_EDGE_LABELED, false);
        check_gnp(2, 3, 0, IGRAPH_DIRECTED, modes[i], IGRAPH_MULTIPLE, IGRAPH_EDGE_LABELED, false);

        check_gnp(8, 15, 0.8, IGRAPH_DIRECTED, modes[i], IGRAPH_NO_MULTIPLE, IGRAPH_EDGE_UNLABELED, true);
        check_gnp(6, 3, 1, IGRAPH_DIRECTED, modes[i], IGRAPH_NO_MULTIPLE, IGRAPH_EDGE_UNLABELED, true);
        check_gnp(8, 15, 0.8, IGRAPH_DIRECTED, modes[i], IGRAPH_MULTIPLE, IGRAPH_EDGE_UNLABELED, true);
        check_gnp(6, 3, 1, IGRAPH_DIRECTED, modes[i], IGRAPH_MULTIPLE, IGRAPH_EDGE_UNLABELED, true);
        check_gnp(8, 15, 0.8, IGRAPH_DIRECTED, modes[i], IGRAPH_MULTIPLE, IGRAPH_EDGE_LABELED, true);
        check_gnp(6, 3, 1, IGRAPH_DIRECTED, modes[i], IGRAPH_MULTIPLE, IGRAPH_EDGE_LABELED, true);
    }

    VERIFY_FINALLY_STACK();

    CHECK_ERROR(igraph_bipartite_game_gnm(&graph, NULL, 0, 10, 20, IGRAPH_DIRECTED, IGRAPH_ALL, IGRAPH_SIMPLE_SW,
                                          false), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_bipartite_game_gnm(&graph, NULL, 0, 10, 20, IGRAPH_DIRECTED, IGRAPH_ALL, IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED), IGRAPH_EINVAL);

    CHECK_ERROR(igraph_bipartite_game_gnm(&graph, NULL, 10, 10, 201, IGRAPH_DIRECTED, IGRAPH_ALL, IGRAPH_SIMPLE_SW,
                                          false), IGRAPH_EINVAL);

    CHECK_ERROR(igraph_bipartite_game_gnm(&graph, NULL, -1, 10, 20, IGRAPH_DIRECTED, IGRAPH_ALL, IGRAPH_SIMPLE_SW,
                                          false), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_bipartite_game_gnm(&graph, NULL, 10, -1, 20, IGRAPH_DIRECTED, IGRAPH_ALL, IGRAPH_SIMPLE_SW,
                                          false), IGRAPH_EINVAL);

    CHECK_ERROR(igraph_bipartite_iea_game(&graph, NULL, 0, 10, 20, IGRAPH_DIRECTED, IGRAPH_ALL), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_bipartite_iea_game(&graph, NULL, 10, 0, 20, IGRAPH_UNDIRECTED, IGRAPH_ALL), IGRAPH_EINVAL);

    CHECK_ERROR(igraph_bipartite_iea_game(&graph, NULL, -1, 10, 20, IGRAPH_DIRECTED, IGRAPH_ALL), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_bipartite_iea_game(&graph, NULL, 10, -1, 20, IGRAPH_UNDIRECTED, IGRAPH_ALL), IGRAPH_EINVAL);

    CHECK_ERROR(igraph_bipartite_game_gnp(&graph, NULL, -1, 10, 0.1, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_SIMPLE_SW,
                                          false), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_bipartite_game_gnp(&graph, NULL, 10, -1, 0.9, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_SIMPLE_SW,
                                          false), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_bipartite_game_gnp(&graph, NULL, 10, 10, 1.1, IGRAPH_UNDIRECTED, IGRAPH_ALL, IGRAPH_SIMPLE_SW,
                                          false), IGRAPH_EINVAL);

    VERIFY_FINALLY_STACK();

    return 0;
}
