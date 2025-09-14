/*
   igraph library.
   Copyright (C) 2008-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

int main(void) {
    igraph_t g;
    igraph_int_t edges_array[] = {0, 1, 1, 2, 3, 4, 5, 6, 6, 5, 1, 4, 1, 6, 0, 3 };
    igraph_vector_int_t edges =
        igraph_vector_int_view(edges_array, sizeof(edges_array) / sizeof(edges_array[0]));
    igraph_vector_bool_t types;

    /* Initialize the library. */
    igraph_setup();

    igraph_vector_bool_init(&types, igraph_vector_int_max(&edges) + 1);
    for (igraph_int_t i = 0; i < igraph_vector_bool_size(&types); i++) {
        VECTOR(types)[i] = i % 2;
    }
    igraph_create_bipartite(&g, &types, &edges, IGRAPH_DIRECTED);
    igraph_write_graph_edgelist(&g, stdout);
    igraph_vector_bool_destroy(&types);
    igraph_destroy(&g);

    return 0;
}
