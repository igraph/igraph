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

#include "layout/layout_internal.h"

#include "test_utilities.inc"

int main() {

    /*******************/
    /* Testing the DLA */
    /*******************/
    long int nodes = 10;
    igraph_i_layout_mergegrid_t grid;
    igraph_vector_t x, y, r;
    long int i;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_vector_init(&x, nodes);
    igraph_vector_init(&y, nodes);
    igraph_vector_init(&r, nodes);
    igraph_i_layout_mergegrid_init(&grid, -5, 5, 100, -5, 5, 100);

    /* radius */
    for (i = 0; i < nodes; i++) {
        VECTOR(r)[i] = rand() / (double)RAND_MAX;
    }
    igraph_vector_sort(&r);

    /* place */
    VECTOR(x)[0] = 0;
    VECTOR(y)[0] = 0;
    igraph_i_layout_merge_place_sphere(&grid, 0, 0, VECTOR(r)[nodes - 1], 0);

    for (i = 1; i < nodes; i++) {
        /*     fprintf(stderr, "%li ", i); */
        igraph_i_layout_merge_dla(&grid, i,
                                  igraph_vector_e_ptr(&x, i),
                                  igraph_vector_e_ptr(&y, i),
                                  VECTOR(r)[nodes - i - 1], 0, 0, 4, 7);
        igraph_i_layout_merge_place_sphere(&grid, VECTOR(x)[i], VECTOR(y)[i],
                                           VECTOR(r)[nodes - i - 1], i);
    }

    /*   for (i=0; i<nodes; i++) {  */
    /*     printf("%f %f\n", VECTOR(x)[i], VECTOR(y)[i]); */
    /*   } */

    /*   print_grid(&grid, stdout); */

    igraph_vector_destroy(&x);
    igraph_vector_destroy(&y);
    igraph_vector_destroy(&r);
    igraph_i_layout_mergegrid_destroy(&grid);

    VERIFY_FINALLY_STACK();

    return 0;
}
