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

void call_and_print(igraph_t *graph) {
    igraph_integer_t mut, asym, null;
    IGRAPH_ASSERT(igraph_dyad_census(graph, &mut, &asym, &null) == IGRAPH_SUCCESS);
    printf("Mutual: %" IGRAPH_PRId " ", mut);
    printf("asymmetric: %" IGRAPH_PRId " ", asym);
    printf("null: %" IGRAPH_PRId "\n\n", null);
}


int main() {
    igraph_t g_0, g_1, g_2, g_lm, g_lmu;

    igraph_small(&g_0, 0, IGRAPH_DIRECTED, -1);
    igraph_small(&g_1, 1, IGRAPH_DIRECTED, -1);
    igraph_small(&g_2, 2, IGRAPH_DIRECTED, 0,1, -1);
    igraph_small(&g_lm, 6, IGRAPH_DIRECTED, 0,1, 0,2, 1,1, 1,3, 2,0, 2,0, 2,3, 3,4, 3,4, -1);
    igraph_small(&g_lmu, 6, IGRAPH_UNDIRECTED, 0,1, 0,2, 1,1, 1,3, 2,0, 2,0, 2,3, 3,4, 3,4, -1);

    printf("No vertices:\n");
    call_and_print(&g_0);

    printf("One vertex:\n");
    call_and_print(&g_1);

    printf("Two vertices:\n");
    call_and_print(&g_2);

    printf("Graph with loops and multiple edges:\n");
    call_and_print(&g_lm);

    printf("Same graph, but undirected:\n");
    call_and_print(&g_lmu);

    igraph_destroy(&g_0);
    igraph_destroy(&g_1);
    igraph_destroy(&g_2);
    igraph_destroy(&g_lm);
    igraph_destroy(&g_lmu);

    VERIFY_FINALLY_STACK();
    return 0;
}
