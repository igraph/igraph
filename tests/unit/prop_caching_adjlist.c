/*
   igraph library.
   Copyright (C) 2021-2026  The igraph development team <igraph@igraph.org>

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

int test_single_loop(void) {
    igraph_t graph;
    igraph_adjlist_t adjlist;
    igraph_bool_t has_multiple;

    igraph_small(&graph, 1, IGRAPH_UNDIRECTED, 0, 0, -1);
    igraph_adjlist_init(&graph, &adjlist, IGRAPH_ALL,
                        IGRAPH_LOOPS_TWICE, IGRAPH_NO_MULTIPLE);
    igraph_adjlist_destroy(&adjlist);

    /* The two incidences of a single loop must not be cached as parallel edges. */
    igraph_has_multiple(&graph, &has_multiple);
    IGRAPH_ASSERT(!has_multiple);

    igraph_destroy(&graph);
    return 0;
}

int test_cache_consistency(void) {
    const igraph_int_t edges[][8] = {
        {0},
        {0, 0},
        {0, 0, 0, 0},
        {0, 1, 0, 1},
        {0, 0, 0, 1},
        {0, 0, 0, 1, 0, 1},
        {0, 0, 0, 0, 0, 1},
        {0, 1, 1, 0},
        {0, 1, 1, 1, 1, 2},
        {0, 1, 1, 1, 1, 1, 1, 2}
    };
    const igraph_int_t lengths[] = {0, 2, 4, 4, 4, 6, 6, 4, 6, 8};
    const igraph_bool_t has_loops[] = {false, true, true, false, true, true, true, false, true, true};
    const igraph_bool_t has_parallel[] = {false, false, true, true, false, true, true, false, false, true};
    const igraph_neimode_t modes[] = {IGRAPH_IN, IGRAPH_OUT, IGRAPH_ALL};
    const igraph_loops_t loops[] = {IGRAPH_NO_LOOPS, IGRAPH_LOOPS_ONCE, IGRAPH_LOOPS_TWICE};
    const igraph_bool_t multiple[] = {IGRAPH_NO_MULTIPLE, IGRAPH_MULTIPLE};
    igraph_vector_int_t neighbors;

    igraph_vector_int_init(&neighbors, 0);
    for (int directed = 0; directed < 2; directed++) {
        for (int g = 0; g < 10; g++) {
            igraph_t graph;
            igraph_vector_int_t edge_view = igraph_vector_int_view(edges[g], lengths[g]);
            /* Reciprocal edges are parallel only when the graph is undirected. */
            igraph_bool_t expected_multiple = has_parallel[g] || (!directed && g == 7);
            igraph_create(&graph, &edge_view, 3, directed);

            for (int mode = 0; mode < 3; mode++) {
                for (int loop = 0; loop < 3; loop++) {
                    for (int multi = 0; multi < 2; multi++) {
                        for (int warm = 0; warm < 3; warm++) {
                            igraph_adjlist_t adjlist;
                            igraph_t copy;
                            igraph_bool_t result;

                            igraph_invalidate_cache(&graph);
                            if (warm >= 1) {
                                igraph_has_multiple(&graph, &result);
                                IGRAPH_ASSERT(result == expected_multiple);
                            }
                            if (warm >= 2) {
                                igraph_has_loop(&graph, &result);
                                IGRAPH_ASSERT(result == has_loops[g]);
                            }

                            igraph_adjlist_init(&graph, &adjlist, modes[mode], loops[loop], multiple[multi]);
                            for (igraph_int_t v = 0; v < 3; v++) {
                                igraph_neighbors(&graph, &neighbors, v, modes[mode], loops[loop], multiple[multi]);
                                IGRAPH_ASSERT(igraph_vector_int_all_e(igraph_adjlist_get(&adjlist, v), &neighbors));
                            }
                            igraph_adjlist_destroy(&adjlist);

                            /* Check the graph and a copy before invalidating the cache. */
                            igraph_copy(&copy, &graph);
                            igraph_has_multiple(&copy, &result);
                            IGRAPH_ASSERT(result == expected_multiple);
                            igraph_destroy(&copy);
                            igraph_has_multiple(&graph, &result);
                            IGRAPH_ASSERT(result == expected_multiple);
                            igraph_has_loop(&graph, &result);
                            IGRAPH_ASSERT(result == has_loops[g]);

                            /* This used to abort when the first call poisoned HAS_MULTI. */
                            igraph_adjlist_init(&graph, &adjlist, modes[mode], IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE);
                            igraph_adjlist_destroy(&adjlist);
                            igraph_has_multiple(&graph, &result);
                            IGRAPH_ASSERT(result == expected_multiple);

                            igraph_invalidate_cache(&graph);
                            igraph_has_multiple(&graph, &result);
                            IGRAPH_ASSERT(result == expected_multiple);
                        }
                    }
                }
            }
            igraph_destroy(&graph);
        }
    }
    igraph_vector_int_destroy(&neighbors);
    return 0;
}

int test_caching(void) {
    igraph_t g_simple, g_loop, g_multiloop, g_multi, g_multi_and_loop;
    char *g_desc[] = {"simple", "loop", "multiloop", "multi", "multi and loop"};
    igraph_adjlist_t adjlist;
    igraph_loops_t loops[] = {IGRAPH_NO_LOOPS, IGRAPH_LOOPS_ONCE, IGRAPH_LOOPS_TWICE};
    igraph_bool_t multiple[] = {IGRAPH_NO_MULTIPLE, IGRAPH_MULTIPLE};
    igraph_neimode_t modes[] = {IGRAPH_OUT, IGRAPH_ALL};

    igraph_vector_int_t edge;
    igraph_int_t vloop[] = {0,0};
    igraph_int_t vmult[] = {0,1};

    igraph_full(&g_simple, 5, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_full(&g_loop, 5, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_full(&g_multiloop, 5, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_full(&g_multi, 5, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_full(&g_multi_and_loop, 5, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

    edge = igraph_vector_int_view(vloop, 2);
    igraph_add_edges(&g_loop, &edge, NULL);
    igraph_add_edges(&g_multiloop, &edge, NULL);
    igraph_add_edges(&g_multiloop, &edge, NULL);
    igraph_add_edges(&g_multi_and_loop, &edge, NULL);

    edge = igraph_vector_int_view(vmult, 2);

    igraph_add_edges(&g_multi, &edge, NULL);
    igraph_add_edges(&g_multi_and_loop, &edge, NULL);

    igraph_t *graphs[] = {&g_simple, &g_loop, &g_multiloop, &g_multi, &g_multi_and_loop};

    for (int g = 0; g < 5; g++) {
        for (int loop = 0; loop < 3; loop++) {
            for (int multi = 0; multi < 2; multi++) {
                for (int mode = 0; mode < 2; mode++) {
                    printf("graph: %s, loop: %d multi: %d, mode: %d\n", g_desc[g], loop, multi, mode);

                    igraph_invalidate_cache(graphs[g]);
                    igraph_adjlist_init(graphs[g], &adjlist, modes[mode], loops[loop], multiple[multi]);
                    if (!igraph_i_property_cache_has(graphs[g], IGRAPH_PROP_HAS_LOOP)) {
                        printf("loop not cached\n");
                    } else {
                        printf("loop cached: %d\n", igraph_i_property_cache_get_bool(graphs[g], IGRAPH_PROP_HAS_LOOP));
                    }
                    if (!igraph_i_property_cache_has(graphs[g], IGRAPH_PROP_HAS_MULTI)) {
                        printf("multi not cached\n");
                    } else {
                        printf("multi cached: %d\n", igraph_i_property_cache_get_bool(graphs[g], IGRAPH_PROP_HAS_MULTI));
                    }
                    igraph_adjlist_destroy(&adjlist);
                }
            }
        }
        igraph_destroy(graphs[g]);
    }

    return 0;
}

int main(void) {

    RUN_TEST(test_single_loop);
    RUN_TEST(test_cache_consistency);
    RUN_TEST(test_caching);

    VERIFY_FINALLY_STACK();

    return 0;
}
