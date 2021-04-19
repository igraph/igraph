/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

int main() {
    igraph_t g;
    igraph_integer_t girth;
    igraph_vector_t v;
    igraph_vector_t circle;
    igraph_real_t chord[] = { 0, 50 };

    igraph_ring(&g, 100, IGRAPH_UNDIRECTED, 0, 1);
    igraph_vector_view(&v, chord, sizeof(chord) / sizeof(igraph_real_t));
    igraph_add_edges(&g, &v, 0);
    igraph_girth(&g, &girth, 0);
    if (girth != 51) {
        return 1;
    }

    igraph_destroy(&g);

    /* Special case: null graph */
    igraph_ring(&g, 0, IGRAPH_UNDIRECTED, 0, 1);
    igraph_vector_init(&circle, 1);
    VECTOR(circle)[0] = 2;
    igraph_girth(&g, &girth, &circle);
    if (girth != 0) {
        return 2;
    }
    if (igraph_vector_size(&circle) != 0) {
        return 3;
    }
    igraph_vector_destroy(&circle);
    igraph_destroy(&g);

    return 0;
}
