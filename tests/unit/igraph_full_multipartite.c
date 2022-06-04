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

int main() {
    igraph_t g;
    igraph_vector_int_t vertex_set;
    igraph_vector_int_t types;

    igraph_vector_int_init(&vertex_set, 3);
    igraph_vector_int_init(&types, 0);

    VECTOR(vertex_set)[0] = 2;
    VECTOR(vertex_set)[1] = 3;
    VECTOR(vertex_set)[2] = 3;

    igraph_full_multipartite(&g, &types, &vertex_set, IGRAPH_DIRECTED, IGRAPH_ALL);

    printf("Edge list:\n");
    igraph_write_graph_edgelist(&g, stdout);

    printf("\nVertex type:\n");
    igraph_vector_int_print(&types);
    
    igraph_vector_int_destroy(&vertex_set);
    igraph_vector_int_destroy(&types);
    igraph_destroy(&g);

    return IGRAPH_SUCCESS;
}
