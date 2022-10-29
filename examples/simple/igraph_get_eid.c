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

int main(void) {
    igraph_t g;
    igraph_integer_t eid;
    igraph_vector_int_t hist;
    igraph_integer_t i;

    /* DIRECTED */

    igraph_star(&g, 10, IGRAPH_STAR_OUT, 0);

    igraph_vector_int_init(&hist, 9);

    for (i = 1; i < 10; i++) {
        igraph_get_eid(&g, &eid, 0, i, IGRAPH_DIRECTED, /*error=*/ true);
        VECTOR(hist)[ eid ] = 1;
    }

    igraph_vector_int_print(&hist);

    igraph_vector_int_destroy(&hist);
    igraph_destroy(&g);

    /* UNDIRECTED */

    igraph_star(&g, 10, IGRAPH_STAR_UNDIRECTED, 0);

    igraph_vector_int_init(&hist, 9);

    for (i = 1; i < 10; i++) {
        igraph_get_eid(&g, &eid, 0, i, IGRAPH_UNDIRECTED, /*error=*/ true);
        VECTOR(hist)[ eid ] += 1;
        igraph_get_eid(&g, &eid, i, 0, IGRAPH_DIRECTED, /*error=*/ true);
        VECTOR(hist)[ eid ] += 1;
    }
    igraph_vector_int_print(&hist);

    igraph_vector_int_destroy(&hist);
    igraph_destroy(&g);

    return 0;
}

