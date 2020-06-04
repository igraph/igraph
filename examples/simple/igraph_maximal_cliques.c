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

#define NODES 1000
#define CLIQUE_SIZE 10
#define NO_CLIQUES 10
#define INT(a) (igraph_rng_get_integer(igraph_rng_default(), 0, (a)))

int permutation(igraph_vector_t *vec) {
    int i, r, tmp;
    for (i = 0; i < CLIQUE_SIZE; i++) {
        r = INT(NODES - 1);
        tmp = VECTOR(*vec)[i];
        VECTOR(*vec)[i] = VECTOR(*vec)[r];
        VECTOR(*vec)[r] = tmp;
    }
    return 0;
}

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

void print_and_destroy_cliques(igraph_vector_ptr_t *cliques) {
    int i;
    sort_cliques(cliques);
    for (i = 0; i < igraph_vector_ptr_size(cliques); i++) {
        igraph_vector_t *v = VECTOR(*cliques)[i];
        igraph_vector_print(v);
        igraph_vector_destroy(v);
        igraph_free(v);
    }
}

int main() {

    igraph_t g, g2, cli;
    igraph_vector_t perm;
    igraph_vector_ptr_t cliques;
    igraph_integer_t no;
    int i;

    igraph_rng_seed(igraph_rng_default(), 42);

    /* Create a graph that has a random component, plus a number of
       relatively small cliques */

    igraph_vector_init_seq(&perm, 0, NODES - 1);
    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, NODES, NODES,
                            /*directed=*/ 0, /*loops=*/ 0);
    igraph_full(&cli, CLIQUE_SIZE, /*directed=*/ 0, /*loops=*/ 0);

    for (i = 0; i < NO_CLIQUES; i++) {
        /* Permute vertices of g */
        permutation(&perm);
        igraph_permute_vertices(&g, &g2, &perm);
        igraph_destroy(&g);
        g = g2;

        /* Add a clique */
        igraph_union(&g2, &g, &cli, /*edge_map1=*/ 0, /*edge_map2=*/ 0);
        igraph_destroy(&g);
        g = g2;
    }
    igraph_simplify(&g, /*multiple=*/ 1, /*loop=*/ 0, /*edge_comb=*/ 0);

    igraph_vector_destroy(&perm);
    igraph_destroy(&cli);

    /* Find the maximal cliques */

    igraph_vector_ptr_init(&cliques, 0);
    igraph_maximal_cliques(&g, &cliques, /*min_size=*/ 3,
                           /*max_size=*/ 0 /*no limit*/);
    igraph_maximal_cliques_count(&g, &no, /*min_size=*/ 3,
                                 /*max_size=*/ 0 /*no limit*/);

    if (no != igraph_vector_ptr_size(&cliques)) {
        return 1;
    }

    /* Print and destroy them */

    print_and_destroy_cliques(&cliques);

    /* Clean up */

    igraph_vector_ptr_destroy(&cliques);
    igraph_destroy(&g);

    /* Build a triangle with a loop (thanks to Emmanuel Navarro) */

    igraph_small(&g, 3, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 0, 0, 0, -1);

    /* Find the maximal cliques */

    igraph_vector_ptr_init(&cliques, 0);
    igraph_maximal_cliques(&g, &cliques, /*min_size=*/ 3,
                           /*max_size=*/ 0 /*no limit*/);
    igraph_maximal_cliques_count(&g, &no, /*min_size=*/ 3,
                                 /*max_size=*/ 0 /*no limit*/);

    if (no != igraph_vector_ptr_size(&cliques)) {
        return 2;
    }

    /* Print and destroy them */

    print_and_destroy_cliques(&cliques);

    /* Clean up */

    igraph_vector_ptr_destroy(&cliques);
    igraph_destroy(&g);

    return 0;
}
