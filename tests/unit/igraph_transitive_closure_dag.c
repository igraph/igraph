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

#include "test_utilities.inc"

int main() {

    igraph_t g, g2;
    igraph_vector_t deg;

    igraph_small(&g, 9, IGRAPH_DIRECTED,
                 8, 7, 7, 6, 6, 3, 6, 0, 3, 2, 3, 1, 5, 0, 4, 1,
                 -1);
    igraph_transitive_closure_dag(&g, &g2);

    if (igraph_vcount(&g2) != igraph_vcount(&g)) {
        return 1;
    }
    if (igraph_ecount(&g2) != 19) {
        return 1;
    }

    igraph_vector_init(&deg, 0);
    igraph_degree(&g2, &deg, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS);
    igraph_vector_print(&deg);
    igraph_degree(&g2, &deg, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS);
    igraph_vector_print(&deg);

    igraph_vector_destroy(&deg);
    igraph_destroy(&g2);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
