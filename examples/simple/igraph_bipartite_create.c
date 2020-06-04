/* -*- mode: C -*-  */
/*
   IGraph library.
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

int main() {

    igraph_real_t edges2[] = {0, 1, 1, 2, 3, 4, 5, 6, 6, 5, 1, 4, 1, 6, 0, 3 };
    igraph_real_t edges3[] = {0, 1, 1, 2, 3, 4, 5, 6, 6, 5, 2, 4, 1, 6, 0, 3 };
    igraph_t g;
    igraph_vector_bool_t types;
    igraph_vector_t edges;
    long int i;
    int ret;

    igraph_vector_view(&edges, edges2, sizeof(edges2) / sizeof(igraph_real_t));
    igraph_vector_bool_init(&types, igraph_vector_max(&edges) + 1);
    for (i = 0; i < igraph_vector_bool_size(&types); i++) {
        VECTOR(types)[i] = i % 2;
    }
    igraph_create_bipartite(&g, &types, &edges, /*directed=*/ 1);
    igraph_write_graph_edgelist(&g, stdout);
    igraph_vector_bool_destroy(&types);
    igraph_destroy(&g);

    /* Error handling */
    igraph_set_error_handler(igraph_error_handler_ignore);

    igraph_vector_view(&edges, edges3, sizeof(edges3) / sizeof(igraph_real_t));
    igraph_vector_bool_init(&types, igraph_vector_max(&edges) + 1);
    for (i = 0; i < igraph_vector_bool_size(&types); i++) {
        VECTOR(types)[i] = i % 2;
    }
    ret = igraph_create_bipartite(&g, &types, &edges, /*directed=*/ 1);
    if (ret != IGRAPH_EINVAL) {
        return 1;
    }
    igraph_vector_bool_destroy(&types);
    igraph_destroy(&g);

    return 0;
}




