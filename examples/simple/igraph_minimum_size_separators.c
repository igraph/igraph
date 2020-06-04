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

int print_and_destroy(igraph_vector_ptr_t *ptr) {
    long int i, n = igraph_vector_ptr_size(ptr);
    for (i = 0; i < n; i++) {
        igraph_vector_t *v = VECTOR(*ptr)[i];
        igraph_vector_print(v);
        igraph_vector_destroy(v);
        igraph_free(v);
    }
    igraph_vector_ptr_destroy(ptr);
    return 0;
}

int main() {
    igraph_t g, g2;
    igraph_vector_ptr_t sep;
    igraph_vs_t vs;

    igraph_small(&g, 7, IGRAPH_UNDIRECTED,
                 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0,
                 -1);
    igraph_vector_ptr_init(&sep, 0);
    igraph_minimum_size_separators(&g, &sep);
    print_and_destroy(&sep);
    igraph_destroy(&g);

    /* ----------------------------------------------------------- */

    igraph_small(&g, 5, IGRAPH_UNDIRECTED,
                 0, 3, 1, 3, 2, 3,
                 0, 4, 1, 4, 2, 4,
                 -1);
    igraph_vector_ptr_init(&sep, 0);
    igraph_minimum_size_separators(&g, &sep);
    print_and_destroy(&sep);
    igraph_destroy(&g);

    /* ----------------------------------------------------------- */

    igraph_small(&g, 5, IGRAPH_UNDIRECTED,
                 2, 0, 3, 0, 4, 0,
                 2, 1, 3, 1, 4, 1,
                 -1);
    igraph_vector_ptr_init(&sep, 0);
    igraph_minimum_size_separators(&g, &sep);
    print_and_destroy(&sep);
    igraph_destroy(&g);

    /* ----------------------------------------------------------- */

    igraph_small(&g, 10, IGRAPH_UNDIRECTED,
                 0, 2, 0, 3, 1, 2, 1, 3, 5, 2, 5, 3, 6, 2, 6, 3,
                 7, 2, 7, 3, 8, 2, 8, 3, 9, 2, 9, 3,
                 2, 4, 4, 3,
                 -1);
    igraph_vector_ptr_init(&sep, 0);
    igraph_minimum_size_separators(&g, &sep);
    print_and_destroy(&sep);
    igraph_destroy(&g);

    /* ----------------------------------------------------------- */

    igraph_full(&g, 4, IGRAPH_UNDIRECTED, /*loops=*/ 0);
    igraph_vector_ptr_init(&sep, 0);
    igraph_minimum_size_separators(&g, &sep);
    print_and_destroy(&sep);
    igraph_destroy(&g);

    /* ----------------------------------------------------------- */

    igraph_small(&g, 23, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 0, 4, 0, 5,
                 1, 2, 1, 3, 1, 4, 1, 6,
                 2, 3, 2, 5, 2, 6,
                 3, 4, 3, 5, 3, 6,
                 4, 5, 4, 6, 4, 20,
                 5, 6,
                 6, 7, 6, 10, 6, 13, 6, 18,
                 7, 8, 7, 10, 7, 13,
                 8, 9,
                 9, 11, 9, 12,
                 10, 11, 10, 13,
                 11, 15,
                 12, 15,
                 13, 14,
                 14, 15,
                 16, 17, 16, 18, 16, 19,
                 17, 19, 17, 20,
                 18, 19, 18, 21, 18, 22,
                 19, 20,
                 20, 21, 20, 22,
                 21, 22,
                 -1);

    igraph_vector_ptr_init(&sep, 0);
    igraph_minimum_size_separators(&g, &sep);
    printf("Orig:\n");
    print_and_destroy(&sep);

    igraph_vector_ptr_init(&sep, 0);
    igraph_vs_vector_small(&vs, 0, 1, 2, 3, 4, 5, 6, 16, 17, 18, 19, 20, 21, 22, -1);
    igraph_induced_subgraph(&g, &g2, vs, IGRAPH_SUBGRAPH_AUTO);
    igraph_minimum_size_separators(&g2, &sep);
    printf("1-7,17-23:\n");
    print_and_destroy(&sep);
    igraph_vs_destroy(&vs);
    igraph_destroy(&g2);

    igraph_vector_ptr_init(&sep, 0);
    igraph_vs_vector_small(&vs, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, -1);
    igraph_induced_subgraph(&g, &g2, vs, IGRAPH_SUBGRAPH_AUTO);
    igraph_minimum_size_separators(&g2, &sep);
    printf("7-16:\n");
    print_and_destroy(&sep);
    igraph_vs_destroy(&vs);
    igraph_destroy(&g2);

    igraph_vector_ptr_init(&sep, 0);
    igraph_vs_vector_small(&vs, 16, 17, 18, 19, 20, 21, 22, -1);
    igraph_induced_subgraph(&g, &g2, vs, IGRAPH_SUBGRAPH_AUTO);
    igraph_minimum_size_separators(&g2, &sep);
    printf("17-23:\n");
    print_and_destroy(&sep);
    igraph_vs_destroy(&vs);
    igraph_destroy(&g2);

    igraph_vector_ptr_init(&sep, 0);
    igraph_vs_vector_small(&vs, 6, 7, 10, 13, -1);
    igraph_induced_subgraph(&g, &g2, vs, IGRAPH_SUBGRAPH_AUTO);
    igraph_minimum_size_separators(&g2, &sep);
    printf("7,8,11,14:\n");
    print_and_destroy(&sep);
    igraph_vs_destroy(&vs);
    igraph_destroy(&g2);

    igraph_vector_ptr_init(&sep, 0);
    igraph_vs_vector_small(&vs, 0, 1, 2, 3, 4, 5, 6, -1);
    igraph_induced_subgraph(&g, &g2, vs, IGRAPH_SUBGRAPH_AUTO);
    igraph_minimum_size_separators(&g2, &sep);
    printf("1-7:\n");
    print_and_destroy(&sep);
    igraph_vs_destroy(&vs);
    igraph_destroy(&g2);

    igraph_destroy(&g);

    return 0;
}
