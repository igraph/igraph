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
#include "test_utilities.h"

void print_and_destroy(igraph_t *g) {
    igraph_adjlist_t al;
    igraph_adjlist_init(g, &al, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE);
    printf("Before:\n");
    igraph_adjlist_print(&al);
    igraph_adjlist_simplify(&al);
    printf("After:\n");
    igraph_adjlist_print(&al);
    igraph_destroy(g);
    igraph_adjlist_destroy(&al);
}

int main(void) {
    igraph_t g;

    printf("Graph with no vertices:\n");
    igraph_small(&g, 0, 1, -1);
    print_and_destroy(&g);

    printf("\nGraph with loops and multiple edges:\n");
    igraph_small(&g, 6, 1, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, -1);
    print_and_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
