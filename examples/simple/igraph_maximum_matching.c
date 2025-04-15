/* -*- mode: C -*-  */
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

int main(void) {
    igraph_t graph;
    igraph_vector_int_t matching;
    igraph_integer_t matching_size;

    igraph_small(&graph, 0, IGRAPH_DIRECTED,
                 1,2, 1,3, 2,4, 3,5,
                 4,5, 5,6, 6,8, 6,9,
                 7,8, 8,9, 8,10, 9,10,
                 -1);

    igraph_vector_int_init(&matching,0);

    igraph_maximum_matching(&graph, &matching_size, NULL, &matching, NULL, 0);

    printf("matching size is: %" IGRAPH_PRId "\n", matching_size);
    printf("matching: ");
    igraph_vector_int_print(&matching);

    igraph_vector_int_destroy(&matching);
    igraph_destroy(&graph);

    return 0;
}
