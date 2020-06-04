/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2012  Gabor Csardi <csardi.gabor@gmail.com>
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


/* This test counts motifs using LAD and compares the results with
 * the RANDESU motif finder */
void test_motifs() {
    igraph_t graph;
    igraph_vector_t randesu_counts, lad_counts;
    igraph_vector_t cut_prob;
    int i, n;
    igraph_bool_t equal;
    igraph_integer_t vcount;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_erdos_renyi_game_gnm(&graph, 40, 400, /* directed = */ 1, /* loops = */ 0);
    vcount = igraph_vcount(&graph);

    /* 3-motifs */

    n = 16; /* there are 16 size-3 directed graphs */

    igraph_vector_init(&lad_counts, n);

    for (i = 0; i < n; i++) {
        igraph_t pattern;
        igraph_vector_ptr_t maps;
        igraph_integer_t nAutomorphisms;

        igraph_isoclass_create(&pattern, 3, i, /* directed = */ 1);
        igraph_vector_ptr_init(&maps, 0);

        igraph_subisomorphic_lad(&pattern, &graph, NULL, NULL, NULL, &maps, /* induced = */ 1, 0);

        igraph_count_subisomorphisms_vf2(&pattern, &pattern, NULL, NULL, NULL, NULL, &nAutomorphisms, NULL, NULL, NULL);

        VECTOR(lad_counts)[i] = igraph_vector_ptr_size(&maps) / nAutomorphisms;

        IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&maps, igraph_vector_destroy);
        igraph_vector_ptr_destroy_all(&maps);

        igraph_destroy(&pattern);
    }

    igraph_vector_init(&cut_prob, 3);
    igraph_vector_init(&randesu_counts, 0);
    igraph_motifs_randesu(&graph, &randesu_counts, 3, &cut_prob);

    equal = 1 /* true */;
    for (i = 0; i < n; i++) {
        if (igraph_is_nan(VECTOR(randesu_counts)[i])) {
            continue;
        }
        if (VECTOR(randesu_counts)[i] != VECTOR(lad_counts)[i]) {
            equal = 0;
            break;
        }
    }

    if (! equal) {
        printf("LAD 3-motif count does not agree with RANDESU.\n");
    }

    if (igraph_vector_sum(&lad_counts) != vcount * (vcount - 1) * (vcount - 2) / 6) {
        printf("Total 3-vertex subgraph count is incorrect.\n");
    }

    igraph_vector_destroy(&randesu_counts);
    igraph_vector_destroy(&lad_counts);
    igraph_vector_destroy(&cut_prob);

    /* 4-motifs */

    n = 218; /* there are 218 size-4 directed graphs */

    igraph_vector_init(&lad_counts, n);

    for (i = 0; i < n; i++) {
        igraph_t pattern;
        igraph_vector_ptr_t maps;
        igraph_integer_t nAutomorphisms;

        igraph_isoclass_create(&pattern, 4, i, /* directed = */ 1);
        igraph_vector_ptr_init(&maps, 0);

        igraph_subisomorphic_lad(&pattern, &graph, NULL, NULL, NULL, &maps, /* induced = */ 1, 0);

        igraph_count_subisomorphisms_vf2(&pattern, &pattern, NULL, NULL, NULL, NULL, &nAutomorphisms, NULL, NULL, NULL);

        VECTOR(lad_counts)[i] = igraph_vector_ptr_size(&maps) / nAutomorphisms;

        IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&maps, igraph_vector_destroy);
        igraph_vector_ptr_destroy_all(&maps);

        igraph_destroy(&pattern);
    }

    igraph_vector_init(&cut_prob, 4);
    igraph_vector_init(&randesu_counts, 0);
    igraph_motifs_randesu(&graph, &randesu_counts, 4, &cut_prob);

    equal = 1 /* true */;
    for (i = 0; i < n; i++) {
        if (igraph_is_nan(VECTOR(randesu_counts)[i])) {
            continue;
        }
        if (VECTOR(randesu_counts)[i] != VECTOR(lad_counts)[i]) {
            equal = 0;
            break;
        }
    }

    if (! equal) {
        printf("LAD 4-motif count does not agree with RANDESU.\n");
    }

    if (igraph_vector_sum(&lad_counts) != vcount * (vcount - 1) * (vcount - 2) * (vcount - 3) / 24) {
        printf("Total 4-vertex subgraph count is incorrect.\n");
    }

    igraph_vector_destroy(&randesu_counts);
    igraph_vector_destroy(&lad_counts);
    igraph_vector_destroy(&cut_prob);

    igraph_destroy(&graph);
}


int main() {
    igraph_t target, pattern;
    igraph_bool_t iso;
    igraph_vector_t map;
    igraph_vector_ptr_t maps;
    int i, n, result;
    int domainsvec[] = { 0, 2, 8, -1,
                         4, 5, 6, 7, -1,
                         1, 3, 5, 6, 7, 8, -1,
                         0, 2, 8, -1,
                         1, 3, 7, 8, -1, -2
                       };
    igraph_vector_ptr_t domains;
    igraph_vector_t *v = 0;

    igraph_small(&target, 9, IGRAPH_UNDIRECTED,
                 0, 1, 0, 4, 0, 6,
                 1, 0, 1, 4, 1, 2,
                 2, 1, 2, 3,
                 3, 2, 3, 4, 3, 5, 3, 7, 3, 8,
                 4, 0, 4, 1, 4, 3, 4, 5, 4, 6,
                 5, 6, 5, 4, 5, 3, 5, 8,
                 6, 0, 6, 4, 6, 5,
                 7, 3, 7, 8,
                 8, 5, 8, 3, 8, 7,
                 -1);
    igraph_simplify(&target, /*multiple=*/ 1, /*loops=*/ 0, /*edge_comb=*/ 0);

    igraph_small(&pattern, 5, IGRAPH_UNDIRECTED,
                 0, 1, 0, 4,
                 1, 0, 1, 4, 1, 2,
                 2, 1, 2, 3,
                 3, 2, 3, 4,
                 4, 3, 4, 1, 4, 0,
                 -1);
    igraph_simplify(&pattern, /*multiple=*/ 1, /*loops=*/ 0, /*edge_comb=*/ 0);

    igraph_vector_init(&map, 0);
    igraph_vector_ptr_init(&maps, 0);

    igraph_subisomorphic_lad(&pattern, &target, /*domains=*/ 0, &iso, &map,
                             &maps, /*induced=*/ 0, /*time_limit=*/ 0);

    if (!iso) {
        return 1;
    }
    igraph_vector_print(&map);
    n = igraph_vector_ptr_size(&maps);
    for (i = 0; i < n; i++) {
        igraph_vector_t *v = VECTOR(maps)[i];
        igraph_vector_print(v);
        igraph_vector_destroy(v);
        igraph_free(v);
    }

    printf("---------\n");

    igraph_subisomorphic_lad(&pattern, &target, /*domains=*/ 0, &iso, &map,
                             &maps, /*induced=*/ 1, /*time_limit=*/ 0);

    if (!iso) {
        return 2;
    }
    igraph_vector_print(&map);
    n = igraph_vector_ptr_size(&maps);
    for (i = 0; i < n; i++) {
        igraph_vector_t *v = VECTOR(maps)[i];
        igraph_vector_print(v);
        igraph_vector_destroy(v);
        igraph_free(v);
    }

    printf("---------\n");

    igraph_vector_ptr_init(&domains, 0);
    i = 0;
    while (1) {
        if (domainsvec[i] == -2) {
            break;
        } else if (domainsvec[i] == -1) {
            igraph_vector_ptr_push_back(&domains, v);
            v = 0;
        } else {
            if (!v) {
                v = (igraph_vector_t *) malloc(sizeof(igraph_vector_t));
                igraph_vector_init(v, 0);
            }
            igraph_vector_push_back(v, domainsvec[i]);
        }
        i++;
    }

    igraph_subisomorphic_lad(&pattern, &target, &domains, &iso, &map, &maps,
                             /*induced=*/ 0, /*time_limit=*/ 0);

    if (!iso) {
        return 3;
    }
    igraph_vector_print(&map);
    n = igraph_vector_ptr_size(&maps);
    for (i = 0; i < n; i++) {
        igraph_vector_t *v = VECTOR(maps)[i];
        igraph_vector_print(v);
        igraph_vector_destroy(v);
        igraph_free(v);
    }

    n = igraph_vector_ptr_size(&domains);
    for (i = 0; i < n; i++) {
        igraph_vector_t *v = VECTOR(domains)[i];
        igraph_vector_destroy(v);
        free(v);
    }

    igraph_vector_ptr_destroy(&domains);
    igraph_vector_destroy(&map);
    igraph_vector_ptr_destroy(&maps);

    igraph_destroy(&pattern);
    igraph_destroy(&target);

    printf("---------\n");

    igraph_vector_init(&map, 0);
    igraph_vector_ptr_init(&maps, 0);

    igraph_small(&target, 9, IGRAPH_UNDIRECTED,
                 0, 1, 0, 4, 0, 6,
                 1, 0, 1, 4, 1, 2,
                 2, 1, 2, 3,
                 3, 2, 3, 4, 3, 5, 3, 7, 3, 8,
                 4, 0, 4, 1, 4, 3, 4, 5, 4, 6,
                 5, 6, 5, 4, 5, 3, 5, 8,
                 6, 0, 6, 4, 6, 5,
                 7, 3, 7, 8,
                 8, 5, 8, 3, 8, 7,
                 -1);
    igraph_simplify(&target, /*multiple=*/ 1, /*loops=*/ 0, /*edge_comb=*/ 0);

    igraph_small(&pattern, 0, IGRAPH_DIRECTED, -1);
    igraph_set_error_handler(igraph_error_handler_ignore);
    result = igraph_subisomorphic_lad(&pattern, &target, /*domains=*/ 0,
                                      &iso, &map, &maps, /*induced=*/ 0,
                                      /*time_limit=*/ 0);
    igraph_set_error_handler(igraph_error_handler_abort);
    if (result != IGRAPH_EINVAL) {
        return 4;
    }
    igraph_destroy(&pattern);

    igraph_small(&pattern, 0, IGRAPH_UNDIRECTED, -1);
    igraph_subisomorphic_lad(&pattern, &target, /*domains=*/ 0, &iso, &map, &maps,
                             /*induced=*/ 0, /*time_limit=*/ 0);
    if (!iso) {
        return 5;
    }
    if (igraph_vector_size(&map) != 0) {
        return 6;
    }
    if (igraph_vector_ptr_size(&maps) != 0) {
        return 7;
    }

    igraph_destroy(&pattern);
    igraph_destroy(&target);

    igraph_vector_destroy(&map);
    igraph_vector_ptr_destroy(&maps);

    test_motifs();

    return 0;
}
