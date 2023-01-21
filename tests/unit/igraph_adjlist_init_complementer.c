/*
   IGraph library.
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

    printf("Graph with loops and mutltiple edges, IGRAPH_IN, loops = 0:\n");
    igraph_small(&g, 6, IGRAPH_DIRECTED, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, -1);
    igraph_adjlist_init_complementer(&g, &adjlist, IGRAPH_IN, 0);
    igraph_adjlist_print(&adjlist);
    igraph_adjlist_destroy(&adjlist);

    printf("Same graph, loops = 1:\n");
    igraph_adjlist_init_complementer(&g, &adjlist, IGRAPH_IN, 1);
    igraph_adjlist_print(&adjlist);
    igraph_adjlist_destroy(&adjlist);

    printf("Same graph, IGRAPH_OUT:\n");
    igraph_adjlist_init_complementer(&g, &adjlist, IGRAPH_OUT, 0);
    igraph_adjlist_fprint(&adjlist, stdout); /* to check fprint too */
    igraph_adjlist_destroy(&adjlist);

    igraph_destroy(&g);
    VERIFY_FINALLY_STACK();
    return 0;
}
