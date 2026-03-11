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

    RUN_TEST(test_caching);

    VERIFY_FINALLY_STACK();

    return 0;
}
