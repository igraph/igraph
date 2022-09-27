/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2013  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "test_utilities.h"

void sort_cliques(igraph_vector_int_list_t *cliques) {
    igraph_integer_t i, n = igraph_vector_int_list_size(cliques);
    for (i = 0; i < n; i++) {
        igraph_vector_int_sort(igraph_vector_int_list_get_ptr(cliques, i));
    }
    igraph_vector_int_list_sort(cliques, igraph_vector_int_lex_cmp);
}

void print_and_destroy(igraph_vector_int_list_t *cliques) {
    sort_cliques(cliques);
    print_vector_int_list(cliques);
    igraph_vector_int_list_destroy(cliques);
}

int main(void) {
    igraph_t graph;
    igraph_vector_int_list_t cliques;

    igraph_rng_seed(igraph_rng_default(), 41);
    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNP,
                            /*n=*/ 100, /*p=*/ 0.7, /*directed=*/ 0,
                            /*loops=*/ 0);

    igraph_vector_int_list_init(&cliques, 0);

    igraph_maximal_cliques(&graph, &cliques, /*min_size=*/ 15,
                           /*max_size=*/ 0);

    print_and_destroy(&cliques);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
