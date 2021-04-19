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

int main () {

    igraph_t g;
    igraph_vs_t vs;
    igraph_vit_t vit;
    igraph_integer_t size;

    /* empty graph, all vertices */
    igraph_empty(&g, 10, IGRAPH_DIRECTED);
    igraph_vs_nonadj(&vs, 0, IGRAPH_ALL);
    igraph_vs_size(&g, &vs, &size);
    printf("%li ", (long int) size);
    igraph_vit_create(&g, vs, &vit);
    while (!IGRAPH_VIT_END(vit)) {
        printf("%li ", (long int) IGRAPH_VIT_GET(vit));
        IGRAPH_VIT_NEXT(vit);
    }
    printf("\n");

    igraph_vit_destroy(&vit);
    igraph_vs_destroy(&vs);
    igraph_destroy(&g);

    /* full graph, no vertices */
    igraph_full(&g, 10, IGRAPH_UNDIRECTED, IGRAPH_LOOPS);
    igraph_vs_nonadj(&vs, 0, IGRAPH_ALL);
    igraph_vit_create(&g, vs, &vit);
    while (!IGRAPH_VIT_END(vit)) {
        printf("%li ", (long int) IGRAPH_VIT_GET(vit));
        IGRAPH_VIT_NEXT(vit);
    }
    printf("\n");

    igraph_vit_destroy(&vit);
    igraph_vs_destroy(&vs);
    igraph_destroy(&g);


    return 0;
}
