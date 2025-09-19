/*
   igraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

int main(void) {
    igraph_t g;
    igraph_adjlist_t adjlist;

    printf("Graph with loops and multiple edges, IGRAPH_IN, loops = IGRAPH_NO_LOOPS:\n");
    igraph_small(&g, 6, IGRAPH_DIRECTED, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, -1);
    igraph_adjlist_init_complementer(&g, &adjlist, IGRAPH_IN, IGRAPH_NO_LOOPS);
    igraph_adjlist_print(&adjlist);
    igraph_adjlist_destroy(&adjlist);

    printf("Same graph, IGRAPH_IN, loops = IGRAPH_LOOPS_ONCE:\n");
    igraph_adjlist_init_complementer(&g, &adjlist, IGRAPH_IN, IGRAPH_LOOPS_ONCE);
    igraph_adjlist_print(&adjlist);
    igraph_adjlist_destroy(&adjlist);

    printf("Same graph, IGRAPH_IN, loops = IGRAPH_LOOPS_TWICE (ignored because directed):\n");
    igraph_adjlist_init_complementer(&g, &adjlist, IGRAPH_IN, IGRAPH_LOOPS_TWICE);
    igraph_adjlist_print(&adjlist);
    igraph_adjlist_destroy(&adjlist);

    printf("Same graph, IGRAPH_OUT, loops = IGRAPH_NO_LOOPS:\n");
    igraph_adjlist_init_complementer(&g, &adjlist, IGRAPH_OUT, IGRAPH_NO_LOOPS);
    igraph_adjlist_fprint(&adjlist, stdout); /* to check fprint too */
    igraph_adjlist_destroy(&adjlist);

    printf("Same graph, IGRAPH_OUT, loops = IGRAPH_LOOPS_ONCE:\n");
    igraph_adjlist_init_complementer(&g, &adjlist, IGRAPH_OUT, IGRAPH_LOOPS_ONCE);
    igraph_adjlist_print(&adjlist);
    igraph_adjlist_destroy(&adjlist);

    printf("Same graph, IGRAPH_OUT, loops = IGRAPH_LOOPS_TWICE (ignored because directed):\n");
    igraph_adjlist_init_complementer(&g, &adjlist, IGRAPH_OUT, IGRAPH_LOOPS_TWICE);
    igraph_adjlist_print(&adjlist);
    igraph_adjlist_destroy(&adjlist);

    printf("Same graph, IGRAPH_ALL, loops = IGRAPH_NO_LOOPS:\n");
    igraph_adjlist_init_complementer(&g, &adjlist, IGRAPH_ALL, IGRAPH_NO_LOOPS);
    igraph_adjlist_print(&adjlist);
    igraph_adjlist_destroy(&adjlist);

    printf("Same graph, IGRAPH_ALL, loops = IGRAPH_LOOPS_ONCE:\n");
    igraph_adjlist_init_complementer(&g, &adjlist, IGRAPH_ALL, IGRAPH_LOOPS_ONCE);
    igraph_adjlist_print(&adjlist);
    igraph_adjlist_destroy(&adjlist);

    printf("Same graph, IGRAPH_ALL, loops = IGRAPH_LOOPS_TWICE:\n");
    igraph_adjlist_init_complementer(&g, &adjlist, IGRAPH_ALL, IGRAPH_LOOPS_TWICE);
    igraph_adjlist_print(&adjlist);
    igraph_adjlist_destroy(&adjlist);

    igraph_destroy(&g);
    VERIFY_FINALLY_STACK();
    return 0;
}
