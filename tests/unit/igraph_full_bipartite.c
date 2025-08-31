/*
   igraph library.
   Copyright (C) 2025  The igraph development team

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

void check_full_bipartite(const igraph_t *g, const igraph_vector_bool_t *types,
                          igraph_int_t n1, igraph_int_t n2,
                          igraph_bool_t directed, igraph_neimode_t mode) {
    igraph_int_t vcount = igraph_vcount(g);
    igraph_int_t ecount = igraph_ecount(g);
    igraph_int_t expected_ecount;
    igraph_bool_t bipartite;

    /* Check vertex count */
    IGRAPH_ASSERT(vcount == n1 + n2);

    /* Check that the graph is bipartite */
    igraph_is_bipartite(g, &bipartite, NULL);
    IGRAPH_ASSERT(bipartite);

    /* Check directedness */
    IGRAPH_ASSERT(igraph_is_directed(g) == directed);

    /* Check types vector if provided */
    if (types) {
        IGRAPH_ASSERT(igraph_vector_bool_size(types) == vcount);
        /* First n1 vertices should be false (type 0), rest should be true (type 1) */
        for (igraph_int_t i = 0; i < n1; i++) {
            IGRAPH_ASSERT(VECTOR(*types)[i] == false);
        }
        for (igraph_int_t i = n1; i < vcount; i++) {
            IGRAPH_ASSERT(VECTOR(*types)[i] == true);
        }
    }

    /* Check edge count */
    if (!directed) {
        expected_ecount = n1 * n2;
    } else if (mode == IGRAPH_OUT || mode == IGRAPH_IN) {
        expected_ecount = n1 * n2;
    } else { /* mode == IGRAPH_ALL */
        expected_ecount = 2 * n1 * n2;
    }
    IGRAPH_ASSERT(ecount == expected_ecount);

    /* Check that edges are only between different partitions */
    for (igraph_int_t i = 0; i < ecount; i++) {
        igraph_int_t from = IGRAPH_FROM(g, i);
        igraph_int_t to = IGRAPH_TO(g, i);

        if (types) {
            IGRAPH_ASSERT(VECTOR(*types)[from] != VECTOR(*types)[to]);
        }

        /* For directed graphs with specific modes, check edge directions */
        if (directed) {
            if (mode == IGRAPH_OUT) {
                IGRAPH_ASSERT(from < n1 && to >= n1);
            } else if (mode == IGRAPH_IN) {
                IGRAPH_ASSERT(from >= n1 && to < n1);
            }
            /* For IGRAPH_ALL, both directions are allowed */
        }
    }
}

int main(void) {
    igraph_t graph;
    igraph_vector_bool_t types;
    igraph_int_t n1, n2;

    igraph_vector_bool_init(&types, 0);

    /* Test 1: Small undirected complete bipartite graph */
    printf("Testing small undirected complete bipartite graph...\n");
    n1 = 3; n2 = 4;
    igraph_full_bipartite(&graph, &types, n1, n2, IGRAPH_UNDIRECTED, IGRAPH_ALL);
    check_full_bipartite(&graph, &types, n1, n2, IGRAPH_UNDIRECTED, IGRAPH_ALL);
    igraph_destroy(&graph);

    /* Test 2: Small directed complete bipartite graph with IGRAPH_OUT */
    printf("Testing small directed complete bipartite graph (OUT mode)...\n");
    n1 = 2; n2 = 3;
    igraph_full_bipartite(&graph, &types, n1, n2, IGRAPH_DIRECTED, IGRAPH_OUT);
    check_full_bipartite(&graph, &types, n1, n2, IGRAPH_DIRECTED, IGRAPH_OUT);
    igraph_destroy(&graph);

    /* Test 3: Small directed complete bipartite graph with IGRAPH_IN */
    printf("Testing small directed complete bipartite graph (IN mode)...\n");
    n1 = 2; n2 = 3;
    igraph_full_bipartite(&graph, &types, n1, n2, IGRAPH_DIRECTED, IGRAPH_IN);
    check_full_bipartite(&graph, &types, n1, n2, IGRAPH_DIRECTED, IGRAPH_IN);
    igraph_destroy(&graph);

    /* Test 4: Small directed complete bipartite graph with IGRAPH_ALL */
    printf("Testing small directed complete bipartite graph (ALL mode)...\n");
    n1 = 2; n2 = 2;
    igraph_full_bipartite(&graph, &types, n1, n2, IGRAPH_DIRECTED, IGRAPH_ALL);
    check_full_bipartite(&graph, &types, n1, n2, IGRAPH_DIRECTED, IGRAPH_ALL);
    igraph_destroy(&graph);

    /* Test 5: Empty first partition */
    printf("Testing empty first partition...\n");
    n1 = 0; n2 = 3;
    igraph_full_bipartite(&graph, &types, n1, n2, IGRAPH_UNDIRECTED, IGRAPH_ALL);
    check_full_bipartite(&graph, &types, n1, n2, IGRAPH_UNDIRECTED, IGRAPH_ALL);
    igraph_destroy(&graph);

    /* Test 6: Empty second partition */
    printf("Testing empty second partition...\n");
    n1 = 3; n2 = 0;
    igraph_full_bipartite(&graph, &types, n1, n2, IGRAPH_UNDIRECTED, IGRAPH_ALL);
    check_full_bipartite(&graph, &types, n1, n2, IGRAPH_UNDIRECTED, IGRAPH_ALL);
    igraph_destroy(&graph);

    /* Test 7: Both partitions empty */
    printf("Testing both partitions empty...\n");
    n1 = 0; n2 = 0;
    igraph_full_bipartite(&graph, &types, n1, n2, IGRAPH_UNDIRECTED, IGRAPH_ALL);
    check_full_bipartite(&graph, &types, n1, n2, IGRAPH_UNDIRECTED, IGRAPH_ALL);
    igraph_destroy(&graph);

    /* Test 8: Singleton partitions */
    printf("Testing singleton partitions...\n");
    n1 = 1; n2 = 1;
    igraph_full_bipartite(&graph, &types, n1, n2, IGRAPH_UNDIRECTED, IGRAPH_ALL);
    check_full_bipartite(&graph, &types, n1, n2, IGRAPH_UNDIRECTED, IGRAPH_ALL);
    igraph_destroy(&graph);

    /* Test 9: Test without types vector */
    printf("Testing without types vector...\n");
    n1 = 3; n2 = 2;
    igraph_full_bipartite(&graph, NULL, n1, n2, IGRAPH_UNDIRECTED, IGRAPH_ALL);
    check_full_bipartite(&graph, NULL, n1, n2, IGRAPH_UNDIRECTED, IGRAPH_ALL);
    igraph_destroy(&graph);

    /* Test 10: Larger graph */
    printf("Testing larger graph...\n");
    n1 = 5; n2 = 6;
    igraph_full_bipartite(&graph, &types, n1, n2, IGRAPH_UNDIRECTED, IGRAPH_ALL);
    check_full_bipartite(&graph, &types, n1, n2, IGRAPH_UNDIRECTED, IGRAPH_ALL);
    igraph_destroy(&graph);

    igraph_vector_bool_destroy(&types);

    /* Test error conditions */
    printf("Testing error conditions...\n");
    /* Test negative n1 parameter */
    CHECK_ERROR(igraph_full_bipartite(&graph, NULL, -1, 5, IGRAPH_UNDIRECTED, IGRAPH_ALL), IGRAPH_EINVAL);
    /* Test negative n2 parameter */
    CHECK_ERROR(igraph_full_bipartite(&graph, NULL, 5, -1, IGRAPH_UNDIRECTED, IGRAPH_ALL), IGRAPH_EINVAL);
    /* Test both negative parameters */
    CHECK_ERROR(igraph_full_bipartite(&graph, NULL, -1, -1, IGRAPH_UNDIRECTED, IGRAPH_ALL), IGRAPH_EINVAL);

    VERIFY_FINALLY_STACK();

    return 0;
}
