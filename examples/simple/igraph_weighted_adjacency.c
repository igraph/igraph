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
    igraph_vector_int_t el;
    igraph_integer_t i, j, n;
    char ch = igraph_is_directed(g) ? '>' : '-';

    igraph_vector_int_init(&el, 0);
    igraph_get_edgelist(g, &el, 0);
    n = igraph_ecount(g);

    for (i = 0, j = 0; i < n; i++, j += 2) {
        printf("%" IGRAPH_PRId " --%c %" IGRAPH_PRId ": %" IGRAPH_PRId "\n",
               VECTOR(el)[j], ch, VECTOR(el)[j + 1], (igraph_integer_t)EAN(g, "weight", i));
    }

    igraph_vector_int_destroy(&el);
}

int main() {
    igraph_t g;
    igraph_matrix_t mat;
    igraph_real_t m[4][4] = {
       { 0, 2, 0, 0 },
       { 1, 0, 0, 1 },
       { 2, 0, 1, 0 },
       { 0, 1, 0, 0 }
    };

    /* Set up an attribute handler */
    igraph_set_attribute_table(&igraph_cattribute_table);

    /* Create a matrix view of the 'm' array. Note that C uses column-major
     * ordering, therefore m[0] is the first column, m[1] is the second etc */
    igraph_matrix_view(&mat, m[0], 4, 4);

    /* Create a graph from the weighted adjacency matrix */
    igraph_weighted_adjacency(&g, &mat, IGRAPH_ADJ_DIRECTED, 0, IGRAPH_LOOPS_ONCE);

    /* Print and destroy the graph. The matrix view does not need to be
     * destroyed */
    print(&g);
    igraph_destroy(&g);

    return 0;
}
