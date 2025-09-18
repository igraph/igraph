/*
   igraph library.
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
    igraph_int_t i, n = igraph_vector_int_list_size(cliques);
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
    igraph_vector_int_list_t cliques, cl1, cl2;
    igraph_vector_int_t v1, v2;
    igraph_int_t n, n1, n2;

    igraph_rng_seed(igraph_rng_default(), 42);
    igraph_erdos_renyi_game_gnp(&graph, /*n=*/ 100, /*p=*/ 0.5, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);

    igraph_vector_int_list_init(&cliques, 0);
    igraph_vector_int_list_init(&cl1, 0);
    igraph_vector_int_list_init(&cl2, 0);

    igraph_maximal_cliques_subset(&graph, /*subset=*/ 0,
                                  &cliques, &n, /*outfile=*/ 0,
                                  /*min_size=*/ 9, /*max_size=*/ IGRAPH_UNLIMITED,
                                  IGRAPH_UNLIMITED);

    igraph_vector_int_init_range(&v1,  0, 13);
    igraph_vector_int_init_range(&v2, 13, 100);
    igraph_maximal_cliques_subset(&graph, &v1, &cl1, &n1, /*outfile=*/ 0,
                                  /*min_size=*/ 9, /*max_size=*/ IGRAPH_UNLIMITED,
                                  IGRAPH_UNLIMITED);
    igraph_maximal_cliques_subset(&graph, &v2, &cl2, &n2, /*outfile=*/ 0,
                                  /*min_size=*/ 9, /*max_size=*/ IGRAPH_UNLIMITED,
                                  IGRAPH_UNLIMITED);

    igraph_vector_int_destroy(&v1);
    igraph_vector_int_destroy(&v2);

    if (n1 + n2 != n) {
        return 1;
    }
    if (n1 != igraph_vector_int_list_size(&cl1)) {
        return 2;
    }
    if (n2 != igraph_vector_int_list_size(&cl2)) {
        return 3;
    }

    print_and_destroy(&cliques);
    printf("---\n");
    print_and_destroy(&cl1);
    printf("+\n");
    print_and_destroy(&cl2);

    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
