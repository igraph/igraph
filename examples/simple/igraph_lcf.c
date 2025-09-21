/*
   igraph library.
   Copyright (C) 2007-2025  The igraph development team <igraph@igraph.org>

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

    igraph_t g1, g2;
    igraph_vector_int_t edges;
    igraph_bool_t iso;

    /* Initialize the library. */
    igraph_setup();

    // Heawood graph through LCF notation: [5, -5]^7
    // The number of vertices is normally the number of shifts
    // multiplied by the number of repeats, in this case 2*7 = 14.
    igraph_lcf_small(&g1,
                     /* n */ 14,
                     /* shifts */ 5, -5,
                     /* repeats */ 7,
                     0);

    printf("edges:\n");
    igraph_vector_int_init(&edges, 0);
    igraph_get_edgelist(&g1, &edges, false);
    igraph_vector_int_print(&edges);
    igraph_vector_int_destroy(&edges);

    // Built-in Heawood graph:
    igraph_famous(&g2, "Heawood");
    igraph_isomorphic(&g1, &g2, &iso);
    printf("isomorphic: %s\n", iso ? "true" : "false");
    igraph_destroy(&g2);

    igraph_destroy(&g1);

    return 0;
}
