/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

int print_mincut(const igraph_t *graph, igraph_real_t value,
                 const igraph_vector_t *partition,
                 const igraph_vector_t *partition2,
                 const igraph_vector_t *cut,
                 const igraph_vector_t *capacity) {

    long int i, nc = igraph_vector_size(cut);
    igraph_bool_t directed = igraph_is_directed(graph);

    printf("mincut value: %g\n", (double) value);
    printf("first partition:  ");
    igraph_vector_print(partition);
    printf("second partition: ");
    igraph_vector_print(partition2);
    printf("edges in the cut: ");
    for (i = 0; i < nc; i++) {
        long int edge = VECTOR(*cut)[i];
        long int from = IGRAPH_FROM(graph, edge);
        long int to  = IGRAPH_TO  (graph, edge);
        if (!directed && from > to) {
            igraph_integer_t tmp = from;
            from = to;
            to = tmp;
        }
        printf("%li-%li (%g), ", from, to, VECTOR(*capacity)[edge]);
    }
    printf("\n");

    return 0;
}

int main() {

    igraph_t g;
    igraph_vector_t weights, partition, partition2, cut;
    igraph_real_t value;

    igraph_vector_init(&partition, 0);
    igraph_vector_init(&partition2, 0);
    igraph_vector_init(&cut, 0);

    /* -------------------------------------------- */

    igraph_small(&g, 0, IGRAPH_UNDIRECTED,
                 0, 1, 0, 4, 1, 2, 1, 4, 1, 5, 2, 3, 2, 6, 3, 6, 3, 7, 4, 5, 5, 6, 6, 7,
                 -1);
    igraph_vector_init_int_end(&weights, -1, 2, 3, 3, 2, 2, 4, 2, 2, 2, 3, 1, 3, -1);

    igraph_mincut(&g, &value, &partition, &partition2, &cut, &weights);
    print_mincut(&g, value, &partition, &partition2, &cut, &weights);

    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    /* -------------------------------------------- */

    igraph_small(&g, 6, IGRAPH_DIRECTED,
                 0, 1, 1, 2, 2, 3, 0, 5, 5, 4, 4, 3, 3, 0, -1);
    igraph_vector_init_int_end(&weights, -1, 3, 1, 2, 10, 1, 3, 2, -1);

    igraph_mincut(&g, &value, &partition, &partition2, &cut, &weights);
    print_mincut(&g, value, &partition, &partition2, &cut, &weights);

    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    /* -------------------------------------------- */

    igraph_small(&g, 5, IGRAPH_DIRECTED,
                 4, 3, 3, 2, 2, 1, 1, 0,
                 -1);
    igraph_vector_init_int_end(&weights, -1, 1, 1, 1, 1, -1);
    igraph_mincut(&g, &value, &partition, &partition2, &cut, &weights);
    print_mincut(&g, value, &partition, &partition2, &cut, &weights);

    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    /* -------------------------------------------- */

    igraph_vector_destroy(&cut);
    igraph_vector_destroy(&partition2);
    igraph_vector_destroy(&partition);

    return 0;
}
