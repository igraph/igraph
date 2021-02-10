/*
   IGraph library.
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

#include "test_utilities.inc"

int test_simple_trees() {
    igraph_t g, g2;
    igraph_adjlist_t adjlist;
    igraph_bool_t iso;

    /* Directed, out */
    igraph_tree(&g, 42, 3, IGRAPH_TREE_OUT);
    igraph_adjlist_init(&g, &adjlist, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);
    igraph_adjlist(&g2, &adjlist, IGRAPH_OUT, /*duplicate=*/ 0);
    igraph_isomorphic(&g, &g2, &iso);
    IGRAPH_ASSERT(iso);
    igraph_adjlist_destroy(&adjlist);
    igraph_destroy(&g2);
    igraph_destroy(&g);

    /* Directed, in */
    igraph_tree(&g, 42, 3, IGRAPH_TREE_OUT);
    igraph_adjlist_init(&g, &adjlist, IGRAPH_IN, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);
    igraph_adjlist(&g2, &adjlist, IGRAPH_IN, /*duplicate=*/ 0);
    igraph_isomorphic(&g, &g2, &iso);
    IGRAPH_ASSERT(iso);
    igraph_adjlist_destroy(&adjlist);
    igraph_destroy(&g2);
    igraph_destroy(&g);

    /* Undirected */
    igraph_tree(&g, 42, 3, IGRAPH_TREE_UNDIRECTED);
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

int test_loop_elimination_for_undirected_graph() {
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

int test_loop_elimination_for_directed_graph() {
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

int test_multiedge_elimination_for_undirected_graph() {
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

int test_multiedge_elimination_for_directed_graph() {
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

int main() {
    int retval;

    RUN_TEST(test_simple_trees);

    RUN_TEST(test_loop_elimination_for_undirected_graph);
    RUN_TEST(test_loop_elimination_for_directed_graph);
    RUN_TEST(test_multiedge_elimination_for_undirected_graph);
    RUN_TEST(test_multiedge_elimination_for_directed_graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
