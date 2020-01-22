/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

void vector_print(igraph_vector_t *v) {
    long int i;
    for (i = 0; i < igraph_vector_size(v); i++) {
        printf(" %li", (long int) VECTOR(*v)[i]);
    }
    printf("\n");
}

int main() {

    igraph_t g;
    igraph_vector_t vids, layers, parents;

    igraph_ring(&g, 10, IGRAPH_UNDIRECTED, 0, 0);
    igraph_vector_init(&vids, 0);
    igraph_vector_init(&layers, 0);
    igraph_vector_init(&parents, 0);
    igraph_i_bfs(&g, 0, IGRAPH_ALL, &vids, &layers, &parents);
    vector_print(&vids);
    vector_print(&layers);
    vector_print(&parents);
    igraph_destroy(&g);

    igraph_tree(&g, 20, 2, IGRAPH_TREE_UNDIRECTED);
    igraph_i_bfs(&g, 0, IGRAPH_ALL, &vids, &layers, &parents);
    vector_print(&vids);
    vector_print(&layers);
    vector_print(&parents);
    igraph_destroy(&g);

    igraph_vector_destroy(&vids);
    igraph_vector_destroy(&layers);
    igraph_vector_destroy(&parents);

    return 0;
}
