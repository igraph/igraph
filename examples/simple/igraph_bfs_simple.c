/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2021  The igraph development team

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

int main(void) {

    igraph_t g;
    igraph_vector_int_t vids, layers, parents;

    igraph_ring(&g, 10, IGRAPH_UNDIRECTED, 0, 0);

    igraph_vector_int_init(&vids, 0);
    igraph_vector_int_init(&layers, 0);
    igraph_vector_int_init(&parents, 0);

    igraph_bfs_simple(&g, 0, IGRAPH_ALL, &vids, &layers, &parents);

    igraph_vector_int_print(&vids);
    igraph_vector_int_print(&layers);
    igraph_vector_int_print(&parents);

    igraph_destroy(&g);

    igraph_vector_int_destroy(&vids);
    igraph_vector_int_destroy(&layers);
    igraph_vector_int_destroy(&parents);

    return 0;
}
