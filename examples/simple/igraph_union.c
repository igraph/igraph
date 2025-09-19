/*
   igraph library.
   Copyright (C) 2006-2023  The igraph development team <igraph@igraph.org>

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

int main(void) {
    igraph_t left, right, uni;
    igraph_vector_int_t edge_map1, edge_map2;

    /* Initialize the library. */
    igraph_setup();

    igraph_vector_int_init(&edge_map1, 0);
    igraph_vector_int_init(&edge_map2, 0);

    igraph_small(&left, 4, IGRAPH_DIRECTED,
                 0,1, 1,2, 2,2, 2,3, 2,3, 3,2,
                 -1);

    igraph_small(&right, 6, IGRAPH_DIRECTED,
                 0,1, 1,2, 2,2, 2,3, 5,2,
                 -1);

    igraph_union(&uni, &left, &right, &edge_map1, &edge_map2);
    igraph_write_graph_edgelist(&uni, stdout);
    igraph_vector_int_print(&edge_map1);
    igraph_vector_int_print(&edge_map2);

    igraph_destroy(&uni);
    igraph_destroy(&left);
    igraph_destroy(&right);

    igraph_vector_int_destroy(&edge_map2);
    igraph_vector_int_destroy(&edge_map1);

    return 0;
}
