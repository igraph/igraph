/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2020  The igraph development team <igraph@igraph.org>

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
#include <stdio.h>

int main() {
    igraph_t graph;
    igraph_vector_t res;

    /* Test graph taken from http://en.wikipedia.org/wiki/Topological_sorting
     * @ 05.03.2006 */
    igraph_small(&graph, 8, IGRAPH_DIRECTED,
                 0, 3, 0, 4, 1, 3, 2, 4, 2, 7, 3, 5, 3, 6, 3, 7, 4, 6,
                 -1);

    igraph_vector_init(&res, 0);

    /* Sort the vertices in "increasing" order. */
    igraph_topological_sorting(&graph, &res, IGRAPH_OUT);
    igraph_vector_print(&res);
    printf("\n");

    /* Sort the vertices in "decreasing" order. */
    igraph_topological_sorting(&graph, &res, IGRAPH_IN);
    igraph_vector_print(&res);

    /* Destroy data structures when done using them. */
    igraph_destroy(&graph);
    igraph_vector_destroy(&res);

    return 0;
}
