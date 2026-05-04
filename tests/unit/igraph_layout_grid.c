/*
   igraph library.
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
#include <stdlib.h>

#include "test_utilities.h"

int main(void) {
    igraph_t g;
    igraph_matrix_t coords;
    igraph_vector_int_t dimvector;

    igraph_empty(&g, 15, 0);
    igraph_matrix_init(&coords, 0, 0);

    /* Predefined width, 2D */
    igraph_layout_grid(&g, &coords, 5);
    igraph_matrix_print(&coords);
    printf("===\n");

    /* Automatic width, 2D */
    igraph_layout_grid(&g, &coords, -1);
    igraph_matrix_print(&coords);
    printf("===\n");

    /* Predefined width and height, 3D */
    igraph_layout_grid_3d(&g, &coords, 4, 2);
    igraph_matrix_print(&coords);
    printf("=====\n");

    /* Predefined width, 3D */
    igraph_layout_grid_3d(&g, &coords, 4, -1);
    igraph_matrix_print(&coords);
    printf("=====\n");

    /* Predefined height, 3D */
    igraph_layout_grid_3d(&g, &coords, -1, 3);
    igraph_matrix_print(&coords);
    printf("=====\n");

    /* Automatic width and height, 3D */
    igraph_layout_grid_3d(&g, &coords, -1, -1);
    igraph_matrix_print(&coords);
    printf("=====\n");

    /* Create 4-dimensional lattice dimvector for igraph_layout_square testing */
    igraph_vector_int_init(&dimvector, 4);
    VECTOR(dimvector)[0] = 4;
    VECTOR(dimvector)[1] = 3;
    VECTOR(dimvector)[2] = 2;
    VECTOR(dimvector)[3] = 2;

    /* 4-dimensional lattice when graph is initialized */
    igraph_empty(&g, 4*3*2*2, IGRAPH_UNDIRECTED);
    igraph_layout_square(&g, &coords, &dimvector);
    igraph_matrix_print(&coords);
    printf("=======\n");
    
    /* 4-dimensional lattice when graph is NULL */
    igraph_layout_square(NULL, &coords, &dimvector);
    igraph_matrix_print(&coords);    
    printf("=======\n");

    igraph_matrix_destroy(&coords);
    igraph_matrix_destroy(&dimvector);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
