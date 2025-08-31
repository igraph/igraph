/*
   igraph library.
   Copyright (C) 2007-2024  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <igraph.h>
#include <stdio.h>

int main(void) {

    igraph_t g;
    igraph_matrix_int_t merges;
    igraph_vector_int_t membership;

    /* Initialize the library. */
    igraph_setup();

    /* Zachary Karate club */
    igraph_small(&g, 0, IGRAPH_UNDIRECTED,
                 0,  1,  0,  2,  0,  3,  0,  4,  0,  5,
                 0,  6,  0,  7,  0,  8,  0, 10,  0, 11,
                 0, 12,  0, 13,  0, 17,  0, 19,  0, 21,
                 0, 31,  1,  2,  1,  3,  1,  7,  1, 13,
                 1, 17,  1, 19,  1, 21,  1, 30,  2,  3,
                 2,  7,  2,  8,  2,  9,  2, 13,  2, 27,
                 2, 28,  2, 32,  3,  7,  3, 12,  3, 13,
                 4,  6,  4, 10,  5,  6,  5, 10,  5, 16,
                 6, 16,  8, 30,  8, 32,  8, 33,  9, 33,
                 13, 33, 14, 32, 14, 33, 15, 32, 15, 33,
                 18, 32, 18, 33, 19, 33, 20, 32, 20, 33,
                 22, 32, 22, 33, 23, 25, 23, 27, 23, 29,
                 23, 32, 23, 33, 24, 25, 24, 27, 24, 31,
                 25, 31, 26, 29, 26, 33, 27, 33, 28, 31,
                 28, 33, 29, 32, 29, 33, 30, 32, 30, 33,
                 31, 32, 31, 33, 32, 33,
                 -1);

    igraph_matrix_int_init(&merges, 0, 0);
    igraph_vector_int_init(&membership, 0);

    igraph_community_leading_eigenvector(&g, /*weights=*/ NULL,
                                         &merges, &membership,
                                         /*steps=*/ 1,
                                         /*options=*/ NULL,
                                         /*modularity=*/ NULL,
                                         /*start=*/ NULL,
                                         /*eigenvalues=*/ NULL,
                                         /*eigenvectors=*/ NULL,
                                         /*history=*/ NULL,
                                         /*callback=*/ NULL,
                                         /*callback_extra=*/ NULL);

    igraph_matrix_int_print(&merges);
    igraph_vector_int_print(&membership);

    printf("\n");

    /* Make all the steps */
    igraph_community_leading_eigenvector(&g, /*weights=*/ NULL,
                                         &merges, &membership,
                                         /*steps=*/ igraph_vcount(&g),
                                         /*options=*/ NULL,
                                         /*modularity=*/ NULL,
                                         /*start=*/ NULL,
                                         /*eigenvalues=*/ NULL,
                                         /*eigenvectors=*/ NULL,
                                         /*history=*/ NULL,
                                         /*callback=*/ NULL,
                                         /*callback_extra=*/ NULL);

    igraph_matrix_int_print(&merges);
    igraph_vector_int_print(&membership);

    igraph_vector_int_destroy(&membership);
    igraph_matrix_int_destroy(&merges);
    igraph_destroy(&g);

    return 0;
}
