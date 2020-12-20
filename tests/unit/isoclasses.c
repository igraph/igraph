/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge MA, 02139 USA

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

#include "test_utilities.inc"

int main() {

    igraph_t g;
    igraph_vector_t edges;
    igraph_vector_t vids;
    igraph_integer_t class;

    igraph_vector_init_int_end(&edges, -1,
                               0, 1, 1, 3, 1, 4, 1, 6, 3, 1,
                               4, 1, 4, 2, 6, 4, 6, 5, 7, 8,
                               8, 7, 7, 9, 9, 7, 8, 9, 9, 8,
                               -1);
    igraph_create(&g, &edges, 0, IGRAPH_DIRECTED);
    igraph_vector_destroy(&edges);

    igraph_vector_init_int_end(&vids, -1, 1, 4, 6, -1);
    igraph_isoclass_subgraph(&g, &vids, &class);
    printf("class: %i\n", (int)class);
    igraph_vector_destroy(&vids);

    igraph_vector_init_int_end(&vids, -1, 0, 1, 3, -1);
    igraph_isoclass_subgraph(&g, &vids, &class);
    printf("class: %i\n", (int)class);
    igraph_vector_destroy(&vids);

    igraph_vector_init_int_end(&vids, -1, 7, 8, 9, -1);
    igraph_isoclass_subgraph(&g, &vids, &class);
    printf("class: %i\n", (int)class);
    igraph_vector_destroy(&vids);

    igraph_vector_init_int_end(&vids, -1, 0, 2, 5, -1);
    igraph_isoclass_subgraph(&g, &vids, &class);
    printf("class: %i\n", (int)class);
    igraph_vector_destroy(&vids);

    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
