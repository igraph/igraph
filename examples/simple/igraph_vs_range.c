/*
   igraph library.
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

int main(void) {

    igraph_vs_t vs;
    igraph_vit_t vit;
    igraph_t g;
    igraph_int_t size;

    /* Initialize the library. */
    igraph_setup();

    igraph_ring(&g, 10, IGRAPH_UNDIRECTED, 0, 1);
    igraph_vs_range(&vs, 0, 10);
    igraph_vit_create(&g, vs, &vit);
    igraph_vs_size(&g, &vs, &size);
    printf("%" IGRAPH_PRId "", size);

    while (!IGRAPH_VIT_END(vit)) {
        printf(" %" IGRAPH_PRId "", IGRAPH_VIT_GET(vit));
        IGRAPH_VIT_NEXT(vit);
    }
    printf("\n");

    igraph_vit_destroy(&vit);
    igraph_vs_destroy(&vs);
    igraph_destroy(&g);

    return 0;
}
