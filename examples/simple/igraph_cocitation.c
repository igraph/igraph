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

void print_matrix(igraph_matrix_t *m, FILE *f) {
    long int i, j;
    for (i = 0; i < igraph_matrix_nrow(m); i++) {
        for (j = 0; j < igraph_matrix_ncol(m); j++) {
            fprintf(f, " %li", (long int) MATRIX(*m, i, j));
        }
        fprintf(f, "\n");
    }
}

int main() {

    igraph_t g;
    igraph_matrix_t m;

    igraph_small(&g, 0, IGRAPH_DIRECTED,
                 0, 1, 2, 1, 2, 0, 3, 0,
                 -1);

    igraph_matrix_init(&m, 0, 0);
    igraph_bibcoupling(&g, &m, igraph_vss_all());
    print_matrix(&m, stdout);

    igraph_cocitation(&g, &m, igraph_vss_all());
    print_matrix(&m, stdout);

    igraph_matrix_destroy(&m);
    igraph_destroy(&g);

    return 0;
}
