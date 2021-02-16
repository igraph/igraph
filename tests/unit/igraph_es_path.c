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
#include <stdlib.h>

#include "test_utilities.inc"

int main() {

    igraph_t g;
    igraph_es_t es;
    igraph_eit_t eit;
    igraph_integer_t size;

    /* DIRECTED */

    igraph_ring(&g, 10, IGRAPH_DIRECTED, 0, 1);
    igraph_es_path_small(&es, IGRAPH_DIRECTED, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, -1);
    igraph_eit_create(&g, es, &eit);
    igraph_es_size(&g, &es, &size);
    while (!IGRAPH_EIT_END(eit)) {
        long int edge = IGRAPH_EIT_GET(eit);
        igraph_integer_t from, to;
        igraph_edge(&g, edge, &from, &to);
        IGRAPH_EIT_NEXT(eit);
        size--;
    }
    IGRAPH_ASSERT(size == 0);

    igraph_eit_destroy(&eit);
    igraph_es_destroy(&es);
    igraph_destroy(&g);

    /* UNDIRECTED */

    igraph_ring(&g, 10, IGRAPH_UNDIRECTED, 0, 1);
    igraph_es_path_small(&es, IGRAPH_DIRECTED,
                         0, 1, 2, 3, 4, 3, 2, 3, 4, 5, 6, 5, 4, 5, 6, 7, 8, 9, 0, 1, 0, 9, -1);
    igraph_eit_create(&g, es, &eit);
    while (!IGRAPH_EIT_END(eit)) {
        long int edge = IGRAPH_EIT_GET(eit);
        igraph_integer_t from, to;
        igraph_edge(&g, edge, &from, &to);
        IGRAPH_EIT_NEXT(eit);
    }

    igraph_eit_destroy(&eit);
    igraph_es_destroy(&es);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
