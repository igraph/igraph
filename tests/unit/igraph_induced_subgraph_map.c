/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

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
#include <stdio.h>

#include "test_utilities.inc"

int main() {
    igraph_t g, sub;
    igraph_vector_t map, invmap;
    igraph_vector_t keep;
    long int i;

    igraph_small(&g, 9, IGRAPH_DIRECTED, 0, 1, 0, 2, 1, 3, 2, 3,
                 1, 4, 4, 2, 1, 5, 5, 2, 1, 6, 6, 2, 1, 7, 7, 2, 1, 8, 8, 2,
                 -1);
    igraph_vector_init(&map, 0);
    igraph_vector_init(&invmap, 0);
    igraph_vector_init(&keep, igraph_vcount(&g));
    for (i = 0; i < igraph_vector_size(&keep); i++) {
        VECTOR(keep)[i] = i;
    }

    igraph_induced_subgraph_map(&g, &sub,
                                igraph_vss_vector(&keep),
                                IGRAPH_SUBGRAPH_COPY_AND_DELETE,
                                &map, &invmap);

    printf("Map:         ");
    igraph_vector_print(&map);
    printf("Inverse map: ");
    igraph_vector_print(&invmap);
    igraph_write_graph_edgelist(&sub, stdout);

    igraph_vector_destroy(&keep);
    igraph_vector_destroy(&map);
    igraph_vector_destroy(&invmap);

    igraph_destroy(&sub);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
