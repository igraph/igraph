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
#include <stdarg.h>

void print(igraph_t *g) {
    igraph_vector_t el;
    long int i, j, n;
    char ch = igraph_is_directed(g) ? '>' : '-';

    igraph_vector_init(&el, 0);
    igraph_get_edgelist(g, &el, 0);
    n = igraph_ecount(g);

    for (i = 0, j = 0; i < n; i++, j += 2) {
        printf("%ld --%c %ld: %ld\n",
               (long)VECTOR(el)[j], ch, (long)VECTOR(el)[j + 1], (long)EAN(g, "weight", i));
    }
    printf("\n");

    igraph_vector_destroy(&el);
}

int main() {
    igraph_t g;
    igraph_matrix_t mat;
    int m[4][4] = { { 0, 1, 2, 0 }, { 2, 0, 0, 1 }, { 0, 0, 1, 0 }, { 0, 1, 0, 0 } };
    long int i, j;

    igraph_matrix_init(&mat, 4, 4);
    for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) {
            MATRIX(mat, i, j) = m[i][j];
        }
    igraph_set_attribute_table(&igraph_cattribute_table);

    /* [ 0 1 2 0 ]
       [ 2 0 0 1 ]
       [ 0 0 1 0 ]
       [ 0 1 0 0 ] */
    igraph_weighted_adjacency(&g, &mat, IGRAPH_ADJ_DIRECTED, 0, /*loops=*/ 1);
    print(&g);
    igraph_destroy(&g);

    /* [ 0 1 2 0 ]
       [ - 0 0 1 ]
       [ - - 1 0 ]
       [ - - - 0 ] */
    igraph_weighted_adjacency(&g, &mat, IGRAPH_ADJ_UPPER, 0, /*loops=*/ 1);
    print(&g);
    igraph_destroy(&g);

    /* [ 0 - - - ]
       [ 2 0 - - ]
       [ 0 0 1 - ]
       [ 0 1 0 0 ] */
    igraph_weighted_adjacency(&g, &mat, IGRAPH_ADJ_LOWER, 0, /*loops=*/ 1);
    print(&g);
    igraph_destroy(&g);

    /* [ 0 1 0 0 ]
       [ 1 0 0 1 ]
       [ 0 0 1 0 ]
       [ 0 1 0 0 ] */
    igraph_weighted_adjacency(&g, &mat, IGRAPH_ADJ_MIN, 0, /*loops=*/ 1);
    print(&g);
    igraph_destroy(&g);

    /* [ 0 2 2 0 ]
       [ 2 0 0 1 ]
       [ 2 0 1 0 ]
       [ 0 1 0 0 ] */
    igraph_weighted_adjacency(&g, &mat, IGRAPH_ADJ_MAX, 0, /*loops=*/ 1);
    print(&g);
    igraph_destroy(&g);

    /* [ 0 3 2 0 ]
       [ 3 0 0 2 ]
       [ 2 0 1 0 ]
       [ 0 2 0 0 ] */
    igraph_weighted_adjacency(&g, &mat, IGRAPH_ADJ_PLUS, 0, /*loops=*/ 1);
    print(&g);
    igraph_destroy(&g);

    igraph_matrix_destroy(&mat);

    if (IGRAPH_FINALLY_STACK_SIZE() != 0) {
        return 1;
    }

    return 0;
}
