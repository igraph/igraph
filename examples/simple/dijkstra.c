/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2008-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

int print_matrix(const igraph_matrix_t *m) {
    long int nrow = igraph_matrix_nrow(m);
    long int ncol = igraph_matrix_ncol(m);
    long int i, j;
    igraph_real_t val;

    for (i = 0; i < nrow; i++) {
        printf("%li:", i);
        for (j = 0; j < ncol; j++) {
            val = MATRIX(*m, i, j);
            if (igraph_is_inf(val)) {
                if (val < 0) {
                    printf("-inf");
                } else {
                    printf(" inf");
                }
            } else {
                printf(" %3.0f", val);
            }
        }
        printf("\n");
    }
    return 0;
}

int main() {

    igraph_t g;
    igraph_vector_t weights;
    igraph_real_t weights_data[] = { 0, 2, 1, 0, 5, 2, 1, 1, 0, 2, 2, 8, 1, 1, 3, 1, 1, 4, 2, 1 };
    igraph_matrix_t res;

    igraph_small(&g, 10, IGRAPH_DIRECTED,
                 0, 1, 0, 2, 0, 3,    1, 2, 1, 4, 1, 5,
                 2, 3, 2, 6,         3, 2, 3, 6,
                 4, 5, 4, 7,         5, 6, 5, 8, 5, 9,
                 7, 5, 7, 8,         8, 9,
                 5, 2,
                 2, 1,
                 -1);

    igraph_vector_view(&weights, weights_data,
                       sizeof(weights_data) / sizeof(igraph_real_t));

    igraph_matrix_init(&res, 0, 0);
    igraph_shortest_paths_dijkstra(&g, &res, igraph_vss_all(), igraph_vss_all(),
                                   &weights, IGRAPH_OUT);
    print_matrix(&res);

    igraph_matrix_destroy(&res);
    igraph_destroy(&g);

    return 0;
}
