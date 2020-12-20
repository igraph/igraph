/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

#include "test_utilities.inc"

int main() {
    igraph_t graph;
    igraph_matrix_t coords;
    int i;

    igraph_matrix_init(&coords, 0, 0);

    for (i = 0; i < 10; i++) {
        igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNP, /*n=*/ 100,
                                /*p=*/ 2.0 / 100, IGRAPH_UNDIRECTED, /*loops=*/ 0);
        igraph_layout_mds(&graph, &coords, /*dist=*/ 0, /*dim=*/ 2);
        igraph_destroy(&graph);
    }

    igraph_matrix_destroy(&coords);

    VERIFY_FINALLY_STACK();

    return 0;
}
