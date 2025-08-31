/*
   igraph library.
   Copyright (C) 2020-2022  The igraph development team <igraph@igraph.org>

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

#define TEST_GRAPH(name) \
    igraph_count_automorphisms_bliss(&graph, NULL, IGRAPH_BLISS_F, &info); \
    printf("%s: %s\n", name, info.group_size); \
    igraph_free(info.group_size); \
    igraph_destroy(&graph);

#define TEST_FAMOUS(name) \
    igraph_famous(&graph, name); \
    TEST_GRAPH(name);

int main(void) {
    igraph_t graph;
    igraph_bliss_info_t info;

    TEST_FAMOUS("Frucht");
    TEST_FAMOUS("Coxeter");
    TEST_FAMOUS("Petersen");
    TEST_FAMOUS("Meredith");

    igraph_full(&graph, 23, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    TEST_GRAPH("Complete 23");

    igraph_star(&graph, 17, IGRAPH_STAR_OUT, 0);
    TEST_GRAPH("Directed star 17");

    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    TEST_GRAPH("Null graph");

    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    TEST_GRAPH("Singleton graph");

    VERIFY_FINALLY_STACK();

    return 0;
}
