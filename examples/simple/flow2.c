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

int main(void) {

    igraph_t g;
    igraph_real_t flow_value;
    igraph_vector_int_t cut;
    igraph_vector_t capacity;
    igraph_vector_int_t partition, partition2;
    igraph_vector_t flow;
    igraph_integer_t i;
    igraph_maxflow_stats_t stats;

    igraph_small(&g, 6, IGRAPH_DIRECTED,
                 0, 1, 1, 2, 2, 3, 0, 5, 5, 4, 4, 3, 3, 0, -1);
    igraph_vector_init_int_end(&capacity, -1, 3, 1, 2, 10, 1, 3, 2, -1);
    igraph_vector_int_init(&cut, 0);
    igraph_vector_int_init(&partition, 0);
    igraph_vector_int_init(&partition2, 0);
    igraph_vector_init(&flow, 0);

    igraph_maxflow(&g, &flow_value, &flow, &cut, &partition, &partition2,
                   /*source=*/ 0, /*target=*/ 2, &capacity, &stats);

    igraph_integer_t nc = igraph_vector_int_size(&cut);
    printf("flow value: %g\n", (double) flow_value);
    printf("flow: ");
    igraph_vector_print(&flow);
    printf("first partition:  ");
    igraph_vector_int_print(&partition);
    printf("second partition: ");
    igraph_vector_int_print(&partition2);
    printf("edges in the cut: ");
    for (i = 0; i < nc; i++) {
        igraph_integer_t edge = VECTOR(cut)[i];
        igraph_integer_t from = IGRAPH_FROM(&g, edge);
        igraph_integer_t to  = IGRAPH_TO(&g, edge);
        printf("%" IGRAPH_PRId "-%" IGRAPH_PRId " (%g), ", from, to, VECTOR(capacity)[edge]);
    }
    printf("\n");

    igraph_vector_int_destroy(&cut);
    igraph_vector_int_destroy(&partition2);
    igraph_vector_int_destroy(&partition);
    igraph_vector_destroy(&capacity);
    igraph_vector_destroy(&flow);
    igraph_destroy(&g);

    return 0;
}
