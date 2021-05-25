/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

    igraph_t g;
    igraph_vector_t result;
    long i;

    igraph_vector_init(&result, 0);

    igraph_small(&g, 7, 0, 0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3, 3, 4, 4, 5, 4, 6, 5, 6, -1);
    igraph_convergence_degree(&g, &result, 0, 0);
    for (i = 0; i < igraph_ecount(&g); i++) {
        printf("%.4f ", (float)igraph_vector_e(&result, i));
    }
    printf("\n");
    igraph_destroy(&g);

    igraph_small(&g, 6, 1, 1, 0, 2, 0, 3, 0, 4, 0, 0, 5, -1);
    igraph_convergence_degree(&g, &result, 0, 0);
    for (i = 0; i < igraph_ecount(&g); i++) {
        printf("%.4f ", (float)igraph_vector_e(&result, i));
    }
    printf("\n");
    igraph_destroy(&g);

    igraph_vector_destroy(&result);

    VERIFY_FINALLY_STACK();

    return 0;
}
