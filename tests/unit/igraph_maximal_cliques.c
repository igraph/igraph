/*
   igraph library.
   Copyright (C) 2023  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <igraph.h>
#include "test_utilities.h"

#define NODES 1000
#define CLIQUE_SIZE 10
#define NO_CLIQUES 10

void permutation(igraph_vector_int_t *vec) {
    igraph_int_t i, r, tmp;
    for (i = 0; i < CLIQUE_SIZE; i++) {
        r = RNG_INTEGER(0, NODES - 1);
        tmp = VECTOR(*vec)[i];
        VECTOR(*vec)[i] = VECTOR(*vec)[r];
        VECTOR(*vec)[r] = tmp;
    }
}

int sort_cmp(const igraph_vector_int_t *a, const igraph_vector_int_t *b) {
    igraph_int_t i, alen = igraph_vector_int_size(a), blen = igraph_vector_int_size(b);
    if (alen != blen) {
        return (alen < blen) - (alen > blen);
    }
    for (i = 0; i < alen; i++) {
        igraph_int_t ea = VECTOR(*a)[i], eb = VECTOR(*b)[i];
        if (ea != eb) {
            return (ea > eb) - (ea < eb);
        }
    }
    return 0;
}

void sort_cliques(igraph_vector_int_list_t *cliques) {
    igraph_int_t i, n = igraph_vector_int_list_size(cliques);
    for (i = 0; i < n; i++) {
        igraph_vector_int_t *v = igraph_vector_int_list_get_ptr(cliques, i);
        igraph_vector_int_sort(v);
    }
    igraph_vector_int_list_sort(cliques, sort_cmp);
}

void print_cliques(igraph_vector_int_list_t *cliques) {
    igraph_int_t i;
    sort_cliques(cliques);
    for (i = 0; i < igraph_vector_int_list_size(cliques); i++) {
        igraph_vector_int_t *v = igraph_vector_int_list_get_ptr(cliques, i);
        igraph_vector_int_print(v);
    }
}

int main(void) {

    igraph_t g, g2, cli;
    igraph_vector_int_t perm, inv_perm;
    igraph_vector_int_list_t cliques;
    igraph_int_t no;
    igraph_int_t i;

    igraph_rng_seed(igraph_rng_default(), 42);

    /* Create a graph that has a random component, plus a number of
       relatively small cliques */

    igraph_vector_int_init_range(&perm, 0, NODES);
    igraph_vector_int_init(&inv_perm, NODES);
    igraph_erdos_renyi_game_gnm(&g, NODES, NODES, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    igraph_full(&cli, CLIQUE_SIZE, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

    for (i = 0; i < NO_CLIQUES; i++) {
        /* Generate a permutation and then invert it for sake of compatibility
         * with earlier tests when igraph_permute_vertices() took the inverse
         * permutation */
        permutation(&perm);
        igraph_invert_permutation(&perm, &inv_perm);

        /* Permute vertices of g */
        igraph_permute_vertices(&g, &g2, &inv_perm);
        igraph_destroy(&g);
        g = g2;

        /* Add a clique */
        igraph_union(&g2, &g, &cli, /*edge_map1=*/ NULL, /*edge_map2=*/ NULL);
        igraph_destroy(&g);
        g = g2;
    }
    igraph_simplify(&g, /*remove_multiple=*/ true, /*remove_loop=*/ false, /*edge_comb=*/ NULL);

    igraph_vector_int_destroy(&inv_perm);
    igraph_vector_int_destroy(&perm);
    igraph_destroy(&cli);

    /* Find the maximal cliques */

    igraph_vector_int_list_init(&cliques, 0);
    igraph_maximal_cliques(&g, &cliques, /*min_size=*/ 3,
                           /*max_size=*/ IGRAPH_UNLIMITED,
                           IGRAPH_UNLIMITED);
    igraph_maximal_cliques_count(&g, &no, /*min_size=*/ 3,
                                 /*max_size=*/ 0 /*no limit*/);

    if (no != igraph_vector_int_list_size(&cliques)) {
        return 1;
    }

    /* Print them */

    print_cliques(&cliques);

    /* Clean up */

    igraph_vector_int_list_destroy(&cliques);
    igraph_destroy(&g);

    /* Build a triangle with a loop (thanks to Emmanuel Navarro) */

    igraph_small(&g, 3, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 0, 0, 0, -1);

    /* Find the maximal cliques */

    igraph_vector_int_list_init(&cliques, 0);
    igraph_maximal_cliques(&g, &cliques, /*min_size=*/ 3,
                           /*max_size=*/ IGRAPH_UNLIMITED,
                           IGRAPH_UNLIMITED);
    igraph_maximal_cliques_count(&g, &no, /*min_size=*/ 3,
                                 /*max_size=*/ 0 /*no limit*/);

    if (no != igraph_vector_int_list_size(&cliques)) {
        return 2;
    }

    /* Print them */

    print_cliques(&cliques);

    /* Clean up */

    igraph_vector_int_list_destroy(&cliques);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
