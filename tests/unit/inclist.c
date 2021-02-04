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

#define TEST_INCLIST(label, mode, loops) { \
    igraph_inclist_init(&g, &inclist, mode, loops); \
    printf(label ": "); \
    print_inclist(&inclist); \
    printf("\n"); \
    igraph_inclist_destroy(&inclist); \
}

int test_loop_elimination_for_undirected_graph() {
    igraph_t g;
    igraph_inclist_t inclist;

    printf("Testing loop edge elimination in undirected graph\n\n");

    igraph_small(
        &g, 5, /* directed = */ 0,
        /* edge 0 */ 0, 1,
        /* edge 1 */ 0, 3,
        /* edge 2 */ 1, 2,
        /* edge 3 */ 2, 2,
        /* edge 4 */ 2, 3,
        /* edge 5 */ 3, 0,
        /* edge 6 */ 3, 4,
        /* edge 7 */ 4, 0,
        /* edge 8 */ 4, 4,
        /* edge 9 */ 4, 5,
        /* edge 10 */ 4, 6,
        /* edge 11 */ 4, 4,
        /* edge 12 */ 6, 5,
        -1
    );

    /* We are testing IGRAPH_ALL, IGRAPH_IN and IGRAPH_OUT below; it should
     * make no difference */
    TEST_INCLIST("Loops eliminated", IGRAPH_ALL, IGRAPH_NO_LOOPS);
    TEST_INCLIST("Loops listed once", IGRAPH_IN, IGRAPH_LOOPS_ONCE);
    TEST_INCLIST("Loops listed twice", IGRAPH_OUT, IGRAPH_LOOPS_TWICE);

    igraph_destroy(&g);

    printf("============================================================\n\n");

    return 0;
}

int test_loop_elimination_for_directed_graph() {
    igraph_t g;
    igraph_inclist_t inclist;

    printf("Testing loop edge elimination in directed graph\n\n");

    igraph_small(
        &g, 5, /* directed = */ 1,
        /* edge 0 */ 0, 1,
        /* edge 1 */ 0, 3,
        /* edge 2 */ 1, 2,
        /* edge 3 */ 2, 2,
        /* edge 4 */ 2, 3,
        /* edge 5 */ 3, 0,
        /* edge 6 */ 3, 4,
        /* edge 7 */ 4, 0,
        /* edge 8 */ 4, 4,
        /* edge 9 */ 4, 5,
        /* edge 10 */ 4, 6,
        /* edge 11 */ 4, 4,
        /* edge 12 */ 6, 5,
        -1
    );

    TEST_INCLIST("In-edges, loops eliminated", IGRAPH_IN, IGRAPH_NO_LOOPS);
    TEST_INCLIST("In-edges, loops listed once", IGRAPH_IN, IGRAPH_LOOPS_ONCE);
    TEST_INCLIST("In-edges, loops listed once even if IGRAPH_LOOPS_TWICE is given",
        IGRAPH_IN, IGRAPH_LOOPS_TWICE);

    TEST_INCLIST("Out-edges, loops eliminated", IGRAPH_OUT, IGRAPH_NO_LOOPS);
    TEST_INCLIST("Out-edges, loops listed once", IGRAPH_OUT, IGRAPH_LOOPS_ONCE);
    TEST_INCLIST("Out-edges, loops listed once even if IGRAPH_LOOPS_TWICE is given",
        IGRAPH_OUT, IGRAPH_LOOPS_TWICE);

    TEST_INCLIST("In- and out-edges, loops eliminated", IGRAPH_ALL, IGRAPH_NO_LOOPS);
    TEST_INCLIST("In- and out-edges, loops listed once", IGRAPH_ALL, IGRAPH_LOOPS_ONCE);
    TEST_INCLIST("In- and out-edges, loops listed twice", IGRAPH_ALL, IGRAPH_LOOPS_TWICE);

    igraph_destroy(&g);

    printf("============================================================\n\n");

    return 0;
}

int main() {
    int retval;

    RUN_TEST(test_loop_elimination_for_undirected_graph);
    RUN_TEST(test_loop_elimination_for_directed_graph);

    return 0;
}
