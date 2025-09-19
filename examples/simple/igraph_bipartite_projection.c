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
#include <stdlib.h>

int main(void) {
    igraph_t g, p1, p2;
    igraph_vector_bool_t types;
    igraph_vector_int_t mult1, mult2;

    /* Initialize the library. */
    igraph_setup();

    igraph_vector_int_init(&mult1, 0);
    igraph_vector_int_init(&mult2, 0);

    igraph_small(&g, 0, IGRAPH_UNDIRECTED, 0,1, 1,2, 0,3, 3,2, 2,4, 4,5, -1);
    igraph_vector_bool_init_int(&types, 6, 0, 1, 0, 1, 1, 0);

    igraph_bipartite_projection(&g, &types, &p1, &p2, &mult1, &mult2, /*probe1=*/ 1);

    igraph_write_graph_edgelist(&p1, stdout);
    printf("\n");
    igraph_write_graph_edgelist(&p2, stdout);
    printf("\n");

    igraph_vector_int_print(&mult1);
    igraph_vector_int_print(&mult2);

    igraph_vector_int_destroy(&mult1);
    igraph_vector_int_destroy(&mult2);
    igraph_destroy(&p1);
    igraph_destroy(&p2);
    igraph_destroy(&g);
    igraph_vector_bool_destroy(&types);

    return 0;
}
