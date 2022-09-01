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

void print_maps(igraph_vector_int_t *map, igraph_vector_int_list_t *maps) {
    igraph_integer_t n, i;
    igraph_vector_int_print(map);
    n = igraph_vector_int_list_size(maps);
    for (i = 0; i < n; i++) {
        igraph_vector_int_t *v = igraph_vector_int_list_get_ptr(maps, i);
        igraph_vector_int_print(v);
    }
    igraph_vector_int_list_clear(maps);
}

/* This test counts motifs using LAD and compares the results with
 * the RANDESU motif finder */
void test_k_motifs(const igraph_t *graph, const int k, const int class_count, igraph_bool_t directed) {
    igraph_vector_t randesu_counts, lad_counts;
    igraph_vector_t cut_prob;
    igraph_bool_t equal;
    igraph_integer_t i, n;
    igraph_integer_t vcount;
    igraph_real_t expected_count;

    vcount = igraph_vcount(graph);

    n = class_count;

    igraph_vector_init(&lad_counts, n);

    for (i = 0; i < n; i++) {
        igraph_t pattern;
        igraph_vector_int_list_t maps;
        igraph_integer_t nAutomorphisms;

        igraph_isoclass_create(&pattern, k, i, directed);
        igraph_vector_int_list_init(&maps, 0);

        igraph_subisomorphic_lad(&pattern, graph, NULL, NULL, NULL, &maps, /* induced = */ true, 0);

        igraph_count_subisomorphisms_vf2(&pattern, &pattern, NULL, NULL, NULL, NULL, &nAutomorphisms, NULL, NULL, NULL);

        VECTOR(lad_counts)[i] = igraph_vector_int_list_size(&maps) / nAutomorphisms;

        igraph_vector_int_list_destroy(&maps);

        igraph_destroy(&pattern);
    }

    igraph_vector_init(&cut_prob, k);
    igraph_vector_init(&randesu_counts, 0);
    igraph_motifs_randesu(graph, &randesu_counts, k, &cut_prob);

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
        printf("LAD %s %d-motif count does not agree with RANDESU.\n", directed ? "directed" : "undirected", k);
    }

    expected_count = 1;
    for (i = 0; i < k; i++) {
        expected_count *= (vcount - i);
    }
    for (i = 0; i < k; i++) {
        expected_count /= (i + 1);
    }
    if (igraph_vector_sum(&lad_counts) != expected_count) {
        printf("Total %d-vertex %s subgraph count is incorrect.\n", k, directed ? "directed" : "undirected");
    }

    igraph_vector_destroy(&randesu_counts);
    igraph_vector_destroy(&lad_counts);
    igraph_vector_destroy(&cut_prob);
}

void test_motifs() {
    igraph_t graph;
    igraph_integer_t count;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_erdos_renyi_game_gnm(&graph, 30, 400, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);

    igraph_graph_count(3, IGRAPH_DIRECTED, &count);
    test_k_motifs(&graph, 3, count, IGRAPH_DIRECTED);

    igraph_graph_count(4, IGRAPH_DIRECTED, &count);
    test_k_motifs(&graph, 4, count, IGRAPH_DIRECTED);

    igraph_destroy(&graph);
}

void test_motifs_undirected() {
    igraph_t graph;
    igraph_integer_t count;

    igraph_rng_seed(igraph_rng_default(), 137);

    igraph_erdos_renyi_game_gnm(&graph, 18, 100, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

    igraph_graph_count(3, IGRAPH_UNDIRECTED, &count);
    test_k_motifs(&graph, 3, count, IGRAPH_UNDIRECTED);

    igraph_graph_count(4, IGRAPH_UNDIRECTED, &count);
    test_k_motifs(&graph, 4, count, IGRAPH_UNDIRECTED);

    igraph_destroy(&graph);

    /* Use a smaller graph so that the test would not take too long. */
    igraph_erdos_renyi_game_gnm(&graph, 9, 36, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

    igraph_graph_count(5, IGRAPH_UNDIRECTED, &count);
    test_k_motifs(&graph, 5, count, IGRAPH_UNDIRECTED);

    igraph_graph_count(6, IGRAPH_UNDIRECTED, &count);
    test_k_motifs(&graph, 6, count, IGRAPH_UNDIRECTED);

    igraph_destroy(&graph);
}


int main() {
    igraph_t target, pattern;
    igraph_bool_t iso;
    igraph_vector_int_t map;
    igraph_vector_int_list_t maps;
    igraph_integer_t i;
    int domainsvec[] = { 0, 2, 8, -1,
                         4, 5, 6, 7, -1,
                         1, 3, 5, 6, 7, 8, -1,
                         0, 2, 8, -1,
                         1, 3, 7, 8, -1, -2
                       };
    igraph_vector_int_list_t domains;
    igraph_vector_int_t v;

    igraph_small(&target, 9, IGRAPH_UNDIRECTED,
                 0, 1, 0, 4, 0, 6,
                 1, 4, 1, 2,
                 2, 3,
                 3, 4, 3, 5, 3, 7, 3, 8,
                 4, 5, 4, 6,
                 5, 6, 5, 8,
                 7, 8,
                 -1);

    igraph_small(&pattern, 5, IGRAPH_UNDIRECTED,
                 0, 1, 0, 4,
                 1, 0, 1, 4, 1, 2,
                 2, 1, 2, 3,
                 3, 2, 3, 4,
                 4, 3, 4, 1, 4, 0,
                 -1);
    igraph_simplify(&pattern, /*multiple=*/ true, /*loops=*/ false, /*edge_comb=*/ NULL);

    igraph_vector_int_init(&map, 0);
    igraph_vector_int_list_init(&maps, 0);

    igraph_subisomorphic_lad(&pattern, &target, /*domains=*/ NULL, &iso, &map,
                             &maps, /*induced=*/ false, /*time_limit=*/ 0);

    if (!iso) {
        return 1;
    }
    print_maps(&map, &maps);

    printf("---------\n");

    igraph_subisomorphic_lad(&pattern, &target, /*domains=*/ NULL, &iso, &map,
                             &maps, /*induced=*/ true, /*time_limit=*/ 0);

    if (!iso) {
        return 2;
    }
    print_maps(&map, &maps);

    printf("---------\n");

    igraph_vector_int_list_init(&domains, 0);
    i = 0;
    igraph_vector_int_init(&v, 0);
    while (1) {
        if (domainsvec[i] == -2) {
            break;
        } else if (domainsvec[i] == -1) {
            igraph_vector_int_list_push_back_copy(&domains, &v);
            igraph_vector_int_clear(&v);
        } else {
            igraph_vector_int_push_back(&v, domainsvec[i]);
        }
        i++;
    }
    igraph_vector_int_destroy(&v);

    igraph_subisomorphic_lad(&pattern, &target, &domains, &iso, &map, &maps,
                             /*induced=*/ false, /*time_limit=*/ 0);

    if (!iso) {
        return 3;
    }
    print_maps(&map, &maps);

    igraph_vector_int_list_destroy(&domains);
    igraph_vector_int_destroy(&map);
    igraph_vector_int_list_destroy(&maps);

    igraph_destroy(&pattern);
    igraph_destroy(&target);

    printf("---------\n");

    igraph_vector_int_init(&map, 0);
    igraph_vector_int_list_init(&maps, 0);

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
    igraph_simplify(&target, /*multiple=*/ true, /*loops=*/ false, /*edge_comb=*/ NULL);

    igraph_small(&pattern, 0, IGRAPH_UNDIRECTED, -1);
    igraph_subisomorphic_lad(&pattern, &target, /*domains=*/ NULL, &iso, &map, &maps,
                             /*induced=*/ false, /*time_limit=*/ 0);
    if (!iso) {
        return 5;
    }
    if (igraph_vector_int_size(&map) != 0) {
        return 6;
    }
    if (igraph_vector_int_list_size(&maps) != 0) {
        return 7;
    }

    igraph_destroy(&pattern);
    igraph_destroy(&target);

    igraph_vector_int_destroy(&map);
    igraph_vector_int_list_destroy(&maps);

    test_motifs();
    test_motifs_undirected();

    return 0;
}
