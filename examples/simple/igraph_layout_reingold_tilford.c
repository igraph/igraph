/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
#include <math.h>

int main() {

    igraph_t g;
    FILE *f;
    igraph_matrix_t coords;
    /* long int i, n; */

    f = fopen("igraph_layout_reingold_tilford.in", "r");
    igraph_read_graph_edgelist(&g, f, 0, 1);
    igraph_matrix_init(&coords, 0, 0);
    igraph_layout_reingold_tilford(&g, &coords, IGRAPH_IN, 0, 0);

    /*n=igraph_vcount(&g);
    for (i=0; i<n; i++) {
      printf("%6.3f %6.3f\n", MATRIX(coords, i, 0), MATRIX(coords, i, 1));
    }*/

    igraph_matrix_destroy(&coords);
    igraph_destroy(&g);
    return 0;
}
