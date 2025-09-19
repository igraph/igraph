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

int test_simple_trees(void) {
    igraph_t g, g2;
    igraph_adjlist_t adjlist;
    igraph_bool_t iso;

    /* Directed, out */
    igraph_kary_tree(&g, 42, 3, IGRAPH_TREE_OUT);
    igraph_adjlist_init(&g, &adjlist, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);
    igraph_adjlist(&g2, &adjlist, IGRAPH_OUT, /*duplicate=*/ 0);
    igraph_isomorphic(&g, &g2, &iso);
    IGRAPH_ASSERT(iso);
    igraph_adjlist_destroy(&adjlist);
    igraph_destroy(&g2);
    igraph_destroy(&g);

    /* Directed, in */
    igraph_kary_tree(&g, 42, 3, IGRAPH_TREE_OUT);
    igraph_adjlist_init(&g, &adjlist, IGRAPH_IN, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);
    igraph_adjlist(&g2, &adjlist, IGRAPH_IN, /*duplicate=*/ 0);
    igraph_isomorphic(&g, &g2, &iso);
    IGRAPH_ASSERT(iso);
    igraph_adjlist_destroy(&adjlist);
    igraph_destroy(&g2);
    igraph_destroy(&g);

    /* Undirected */
    igraph_kary_tree(&g, 42, 3, IGRAPH_TREE_UNDIRECTED);
    igraph_adjlist_init(&g, &adjlist, IGRAPH_OUT, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE);
    igraph_adjlist(&g2, &adjlist, IGRAPH_ALL, /*duplicate=*/ 1);
    igraph_isomorphic(&g, &g2, &iso);
    IGRAPH_ASSERT(iso);
    igraph_adjlist_destroy(&adjlist);
    igraph_destroy(&g2);
    igraph_destroy(&g);

    return 0;
}

#define TEST_ADJLIST(label, mode, loops, multiple) { \
    igraph_adjlist_init(&g, &adjlist, mode, loops, multiple); \
    printf(label ": "); \
    print_adjlist(&adjlist); \
    printf("\n"); \
    igraph_adjlist_destroy(&adjlist); \
}

#define TEST_LAZY_ADJLIST(label, mode, loops, multiple) { \
    igraph_lazy_adjlist_init(&g, &lazy_adjlist, mode, loops, multiple); \
    printf(label ": "); \
    print_lazy_adjlist(&lazy_adjlist); \
    printf("\n"); \
    igraph_lazy_adjlist_destroy(&lazy_adjlist); \
}

int test_loop_elimination_for_undirected_graph(void) {
    igraph_t g;
    igraph_adjlist_t adjlist;
    igraph_lazy_adjlist_t lazy_adjlist;

    igraph_small(
        &g, 5, /* directed = */ 0,
        0, 1, 0, 3,
        1, 2,
        2, 2, 2, 3,
        3, 0, 3, 4,
        4, 4, 4, 4,
        -1
    );

    printf("Testing loop edge elimination in undirected graph\n\n");

    /* We are testing IGRAPH_ALL, IGRAPH_IN and IGRAPH_OUT below; it should
     * make no difference */
    TEST_ADJLIST("Loops eliminated", IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE);
    TEST_ADJLIST("Loops listed once", IGRAPH_IN, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);
    TEST_ADJLIST("Loops listed twice", IGRAPH_OUT, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE);

    printf("============================================================\n\n");

    printf("Testing lazy loop edge elimination in undirected graph\n\n");

    /* We are testing IGRAPH_ALL, IGRAPH_IN and IGRAPH_OUT below; it should
     * make no difference */
    TEST_LAZY_ADJLIST("Loops eliminated", IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE);
    TEST_LAZY_ADJLIST("Loops listed once", IGRAPH_IN, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);
    TEST_LAZY_ADJLIST("Loops listed twice", IGRAPH_OUT, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE);

    printf("============================================================\n\n");

    igraph_destroy(&g);

    return 0;
}

int test_loop_elimination_for_directed_graph(void) {
    igraph_t g;
    igraph_adjlist_t adjlist;
    igraph_lazy_adjlist_t lazy_adjlist;

    igraph_small(
        &g, 5, /* directed = */ 1,
        0, 1, 0, 3,
        1, 2,
        2, 2, 2, 3,
        3, 0, 3, 4,
        4, 4, 4, 4,
        -1
    );

    printf("Testing loop edge elimination in directed graph\n\n");

    TEST_ADJLIST("In-edges, loops eliminated", IGRAPH_IN, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE);
    TEST_ADJLIST("In-edges, loops listed once", IGRAPH_IN, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);
    TEST_ADJLIST("In-edges, loops listed once even if IGRAPH_LOOPS_TWICE is given",
        IGRAPH_IN, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE);

    TEST_ADJLIST("Out-edges, loops eliminated", IGRAPH_OUT, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE);
    TEST_ADJLIST("Out-edges, loops listed once", IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);
    TEST_ADJLIST("Out-edges, loops listed once even if IGRAPH_LOOPS_TWICE is given",
        IGRAPH_OUT, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE);

    TEST_ADJLIST("In- and out-edges, loops eliminated", IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE);
    TEST_ADJLIST("In- and out-edges, loops listed once", IGRAPH_ALL, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);
    TEST_ADJLIST("In- and out-edges, loops listed twice", IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE);

    printf("============================================================\n\n");

    printf("Testing lazy loop edge elimination in directed graph\n\n");

    TEST_LAZY_ADJLIST("In-edges, loops eliminated", IGRAPH_IN, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE);
    TEST_LAZY_ADJLIST("In-edges, loops listed once", IGRAPH_IN, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);
    TEST_LAZY_ADJLIST("In-edges, loops listed once even if IGRAPH_LOOPS_TWICE is given",
        IGRAPH_IN, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE);

    TEST_LAZY_ADJLIST("Out-edges, loops eliminated", IGRAPH_OUT, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE);
    TEST_LAZY_ADJLIST("Out-edges, loops listed once", IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);
    TEST_LAZY_ADJLIST("Out-edges, loops listed once even if IGRAPH_LOOPS_TWICE is given",
        IGRAPH_OUT, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE);

    TEST_LAZY_ADJLIST("In- and out-edges, loops eliminated", IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE);
    TEST_LAZY_ADJLIST("In- and out-edges, loops listed once", IGRAPH_ALL, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);
    TEST_LAZY_ADJLIST("In- and out-edges, loops listed twice", IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE);

    printf("============================================================\n\n");

    igraph_destroy(&g);

    return 0;
}

int test_multiedge_elimination_for_undirected_graph(void) {
    igraph_t g;
    igraph_adjlist_t adjlist;
    igraph_lazy_adjlist_t lazy_adjlist;

    igraph_small(
        &g, 5, /* directed = */ 0,
        0, 1, 0, 3, 0, 8,
        1, 2,
        2, 2, 2, 3,
        3, 0, 3, 4,
        4, 4, 4, 4, 4, 5, 4, 5, 4, 5,
        5, 6,
        6, 7, 6, 8,
        8, 0,
        -1
    );

    printf("Testing multiple edge elimination in undirected graph\n\n");

    /* We are testing IGRAPH_ALL, IGRAPH_IN and IGRAPH_OUT below; it should
     * make no difference */
    TEST_ADJLIST("Loops also eliminated", IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE);
    TEST_ADJLIST("Loops listed once", IGRAPH_IN, IGRAPH_LOOPS_ONCE, IGRAPH_NO_MULTIPLE);
    TEST_ADJLIST("Loops listed twice", IGRAPH_OUT, IGRAPH_LOOPS_TWICE, IGRAPH_NO_MULTIPLE);

    printf("============================================================\n\n");

    printf("Testing lazy multiple edge elimination in undirected graph\n\n");

    /* We are testing IGRAPH_ALL, IGRAPH_IN and IGRAPH_OUT below; it should
     * make no difference */
    TEST_LAZY_ADJLIST("Loops also eliminated", IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE);
    TEST_LAZY_ADJLIST("Loops listed once", IGRAPH_IN, IGRAPH_LOOPS_ONCE, IGRAPH_NO_MULTIPLE);
    TEST_LAZY_ADJLIST("Loops listed twice", IGRAPH_OUT, IGRAPH_LOOPS_TWICE, IGRAPH_NO_MULTIPLE);

    printf("============================================================\n\n");

    igraph_destroy(&g);

    return 0;
}

int test_multiedge_elimination_for_directed_graph(void) {
    igraph_t g;
    igraph_adjlist_t adjlist;
    igraph_lazy_adjlist_t lazy_adjlist;

    igraph_small(
        &g, 5, /* directed = */ 1,
        0, 1, 0, 3, 0, 8,
        1, 2,
        2, 2, 2, 3,
        3, 0, 3, 4,
        4, 4, 4, 4, 4, 5, 4, 5, 4, 5,
        5, 6,
        6, 7, 6, 8,
        8, 0,
        -1
    );

    printf("Testing multiple edge elimination in directed graph\n\n");

    TEST_ADJLIST("In-edges, loops also eliminated", IGRAPH_IN, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE);
    TEST_ADJLIST("In-edges, loops listed once", IGRAPH_IN, IGRAPH_LOOPS_ONCE, IGRAPH_NO_MULTIPLE);
    TEST_ADJLIST("In-edges, loops listed once even if IGRAPH_LOOPS_TWICE is given",
        IGRAPH_IN, IGRAPH_LOOPS_TWICE, IGRAPH_NO_MULTIPLE);

    TEST_ADJLIST("Out-edges, loops also eliminated", IGRAPH_OUT, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE);
    TEST_ADJLIST("Out-edges, loops listed once", IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_NO_MULTIPLE);
    TEST_ADJLIST("Out-edges, loops listed once even if IGRAPH_LOOPS_TWICE is given",
        IGRAPH_OUT, IGRAPH_LOOPS_TWICE, IGRAPH_NO_MULTIPLE);

    TEST_ADJLIST("In- and out-edges, loops also eliminated", IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE);
    TEST_ADJLIST("In- and out-edges, loops listed once", IGRAPH_ALL, IGRAPH_LOOPS_ONCE, IGRAPH_NO_MULTIPLE);
    TEST_ADJLIST("In- and out-edges, loops listed twice", IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_NO_MULTIPLE);

    printf("============================================================\n\n");

    printf("Testing lazy multiple edge elimination in directed graph\n\n");

    TEST_LAZY_ADJLIST("In-edges, loops also eliminated", IGRAPH_IN, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE);
    TEST_LAZY_ADJLIST("In-edges, loops listed once", IGRAPH_IN, IGRAPH_LOOPS_ONCE, IGRAPH_NO_MULTIPLE);
    TEST_LAZY_ADJLIST("In-edges, loops listed once even if IGRAPH_LOOPS_TWICE is given",
        IGRAPH_IN, IGRAPH_LOOPS_TWICE, IGRAPH_NO_MULTIPLE);

    TEST_LAZY_ADJLIST("Out-edges, loops also eliminated", IGRAPH_OUT, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE);
    TEST_LAZY_ADJLIST("Out-edges, loops listed once", IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_NO_MULTIPLE);
    TEST_LAZY_ADJLIST("Out-edges, loops listed once even if IGRAPH_LOOPS_TWICE is given",
        IGRAPH_OUT, IGRAPH_LOOPS_TWICE, IGRAPH_NO_MULTIPLE);

    TEST_LAZY_ADJLIST("In- and out-edges, loops also eliminated", IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE);
    TEST_LAZY_ADJLIST("In- and out-edges, loops listed once", IGRAPH_ALL, IGRAPH_LOOPS_ONCE, IGRAPH_NO_MULTIPLE);
    TEST_LAZY_ADJLIST("In- and out-edges, loops listed twice", IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_NO_MULTIPLE);

    printf("============================================================\n\n");

    igraph_destroy(&g);

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

    igraph_full(&g_simple, 5, IGRAPH_UNDIRECTED, /*loops*/ 0);
    igraph_full(&g_loop, 5, IGRAPH_UNDIRECTED, /*loops*/ 0);
    igraph_full(&g_multiloop, 5, IGRAPH_UNDIRECTED, /*loops*/ 0);
    igraph_full(&g_multi, 5, IGRAPH_UNDIRECTED, /*loops*/ 0);
    igraph_full(&g_multi_and_loop, 5, IGRAPH_UNDIRECTED, /*loops*/ 0);

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

    RUN_TEST(test_simple_trees);

    RUN_TEST(test_loop_elimination_for_undirected_graph);
    RUN_TEST(test_loop_elimination_for_directed_graph);
    RUN_TEST(test_multiedge_elimination_for_undirected_graph);
    RUN_TEST(test_multiedge_elimination_for_directed_graph);

    RUN_TEST(test_caching);

    VERIFY_FINALLY_STACK();

    return 0;
}
