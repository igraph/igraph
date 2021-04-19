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

int main() {

    igraph_t g;
    igraph_vector_t v = IGRAPH_VECTOR_NULL;
    igraph_real_t edges[] = { 0, 1, 1, 2, 2, 2, 2, 3, 2, 4, 3, 4 };
    igraph_vector_t v2;
    long int i;
    igraph_vit_t vit;
    igraph_vs_t vs;
    igraph_integer_t size;

    igraph_vector_view(&v, edges, sizeof(edges) / sizeof(igraph_real_t));
    igraph_create(&g, &v, 0, IGRAPH_DIRECTED);

    /* Create iterator based on a vector (view) */
    igraph_vector_init(&v2, 6);
    VECTOR(v2)[0] = 0;
    VECTOR(v2)[1] = 2;
    VECTOR(v2)[2] = 4;
    VECTOR(v2)[3] = 0;
    VECTOR(v2)[4] = 2;
    VECTOR(v2)[5] = 4;

    igraph_vit_create(&g, igraph_vss_vector(&v2), &vit);

    i = 0;
    while (!IGRAPH_VIT_END(vit)) {
        if (IGRAPH_VIT_GET(vit) != VECTOR(v2)[i]) {
            return 1;
        }
        IGRAPH_VIT_NEXT(vit);
        i++;
    }
    if (i != igraph_vector_size(&v2)) {
        return 2;
    }

    igraph_vit_destroy(&vit);
    igraph_vector_destroy(&v2);

    /* Create small vector iterator */

    igraph_vs_vector_small(&vs, 0, 2, 4, 0, 2, 4, 2, -1);
    igraph_vit_create(&g, vs, &vit);
    igraph_vs_size(&g, &vs, &size);
    printf("%li ", (long int) size);
    for (; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
        printf("%li ", (long int) IGRAPH_VIT_GET(vit));
    }
    printf("\n");

    igraph_vit_destroy(&vit);
    igraph_vs_destroy(&vs);

    /* Clean up */

    igraph_destroy(&g);

    return 0;
}
