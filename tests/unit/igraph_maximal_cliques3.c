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

#include "test_utilities.inc"

int sort_cmp(const void *a, const void *b) {
    const igraph_vector_t **da = (const igraph_vector_t **) a;
    const igraph_vector_t **db = (const igraph_vector_t **) b;
    int i, alen = igraph_vector_size(*da), blen = igraph_vector_size(*db);
    if (alen != blen) {
        return (alen < blen) - (alen > blen);
    }
    for (i = 0; i < alen; i++) {
        int ea = VECTOR(**da)[i], eb = VECTOR(**db)[i];
        if (ea != eb) {
            return (ea > eb) - (ea < eb);
        }
    }
    return 0;
}

void sort_cliques(igraph_vector_ptr_t *cliques) {
    int i, n = igraph_vector_ptr_size(cliques);
    for (i = 0; i < n; i++) {
        igraph_vector_t *v = VECTOR(*cliques)[i];
        igraph_vector_sort(v);
    }
    igraph_qsort(VECTOR(*cliques), (size_t) n,
                 sizeof(igraph_vector_t *), sort_cmp);
}

int print_and_destroy(igraph_vector_ptr_t *cliques) {
    int i, n = igraph_vector_ptr_size(cliques);
    sort_cliques(cliques);
    for (i = 0; i < n; i++) {
        igraph_vector_t *v = VECTOR(*cliques)[i];
        igraph_vector_print(v);
        igraph_vector_destroy(v);
    }
    igraph_vector_ptr_destroy_all(cliques);
    return 0;
}

int main() {
    igraph_t graph;
    igraph_vector_ptr_t cliques;

    igraph_rng_seed(igraph_rng_default(), 42);
    igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNP,
                            /*n=*/ 100, /*p=*/ 0.7, /*directed=*/ 0,
                            /*loops=*/ 0);

    igraph_vector_ptr_init(&cliques, 0);

    igraph_maximal_cliques(&graph, &cliques, /*min_size=*/ 15,
                           /*max_size=*/ 0);

    print_and_destroy(&cliques);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
