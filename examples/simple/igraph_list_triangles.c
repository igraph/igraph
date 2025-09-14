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

int main(void) {
    igraph_t g;
    igraph_vector_int_t v;

    /* Initialize the library. */
    igraph_setup();

    igraph_full(&g, 5, 0, IGRAPH_NO_LOOPS);

    printf("Triangles in a full graph of 5 vertices:\n");
    igraph_vector_int_init(&v, 0);
    igraph_list_triangles(&g, &v);

    const igraph_matrix_int_t result = igraph_matrix_int_view_from_vector(&v, /* nrow = */ 3);
    igraph_matrix_int_print(&result);

    igraph_vector_int_destroy(&v);
    igraph_destroy(&g);

    return 0;
}
