/* igraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

/* This test counts motifs using LAD and compares the results with
 * the RANDESU motif finder */
void test_k_motifs(const igraph_t *graph, const int k, const int class_count, igraph_bool_t directed) {
    igraph_vector_t randesu_counts, lad_counts;
    igraph_bool_t equal;
    igraph_int_t i, n;
    igraph_int_t vcount;
    igraph_real_t expected_count;

    vcount = igraph_vcount(graph);

    n = class_count;

    igraph_vector_init(&lad_counts, n);

    for (i = 0; i < n; i++) {
        igraph_t pattern;
        igraph_vector_int_list_t maps;
        igraph_int_t nAutomorphisms;

        igraph_isoclass_create(&pattern, k, i, directed);
        igraph_vector_int_list_init(&maps, 0);

        igraph_subisomorphic_lad(&pattern, graph, NULL, NULL, NULL, &maps, /* induced = */ true);

        igraph_count_subisomorphisms_vf2(&pattern, &pattern, NULL, NULL, NULL, NULL, &nAutomorphisms, NULL, NULL, NULL);

        VECTOR(lad_counts)[i] = igraph_vector_int_list_size(&maps) / nAutomorphisms;

        igraph_vector_int_list_destroy(&maps);

        igraph_destroy(&pattern);
    }

    igraph_vector_init(&randesu_counts, 0);
    igraph_motifs_randesu(graph, &randesu_counts, k, NULL);

    equal = 1 /* true */;
    for (i = 0; i < n; i++) {
        if (isnan(VECTOR(randesu_counts)[i])) {
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
}

void test_motifs(void) {
    igraph_t graph;
    igraph_int_t count;

    igraph_rng_seed(igraph_rng_default(), 42);

    /* The graph is chosen to have approximately 50% density
     * so that most motifs have a high chance of appearing. */
    igraph_erdos_renyi_game_gnm(&graph, 30, 400, IGRAPH_DIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);

    igraph_graph_count(3, IGRAPH_DIRECTED, &count);
    test_k_motifs(&graph, 3, count, IGRAPH_DIRECTED);

    igraph_graph_count(4, IGRAPH_DIRECTED, &count);
    test_k_motifs(&graph, 4, count, IGRAPH_DIRECTED);

    igraph_destroy(&graph);
}

void test_motifs_undirected(void) {
    igraph_t graph;
    igraph_int_t count;

    igraph_rng_seed(igraph_rng_default(), 137);

    /* The graph is chosen to have slightly higher than 50% density
     * so that most connected motifs have a high chance of appearing. */
    igraph_erdos_renyi_game_gnm(&graph, 18, 80, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);

    igraph_graph_count(3, IGRAPH_UNDIRECTED, &count);
    test_k_motifs(&graph, 3, count, IGRAPH_UNDIRECTED);

    igraph_graph_count(4, IGRAPH_UNDIRECTED, &count);
    test_k_motifs(&graph, 4, count, IGRAPH_UNDIRECTED);

    igraph_destroy(&graph);

    /* Use a smaller graph so that the test would not take too long.
     * The graph is chosen to have slightly higher than 50% density
     * so that most connected motifs have a high chance of appearing. */
    igraph_erdos_renyi_game_gnm(&graph, 12, 36, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);

    igraph_graph_count(5, IGRAPH_UNDIRECTED, &count);
    test_k_motifs(&graph, 5, count, IGRAPH_UNDIRECTED);

    igraph_graph_count(6, IGRAPH_UNDIRECTED, &count);
    test_k_motifs(&graph, 6, count, IGRAPH_UNDIRECTED);

    igraph_destroy(&graph);
}


int main(void) {
    igraph_t pattern, target;
    igraph_bool_t iso;
    igraph_vector_int_t map;
    igraph_vector_int_list_t maps;

    igraph_vector_int_init(&map, 0);
    igraph_vector_int_list_init(&maps, 0);

    igraph_small(&target, 9, IGRAPH_UNDIRECTED,
                 0, 1, 0, 4, 0, 6,
                 1, 4, 1, 2,
                 2, 3,
                 3, 4, 3, 5, 3, 7, 3, 8,
                 4, 5, 4, 6,
                 5, 6, 5, 8,
                 7, 8,
                 -1);

    igraph_small(&pattern, 0, IGRAPH_UNDIRECTED, -1);
    igraph_subisomorphic_lad(&pattern, &target, /*domains=*/ NULL, &iso, &map, &maps,
                             /*induced=*/ false);

    IGRAPH_ASSERT(iso);
    IGRAPH_ASSERT(igraph_vector_int_size(&map) == 0);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&maps) == 1);
    IGRAPH_ASSERT(igraph_vector_int_size(igraph_vector_int_list_get_ptr(&maps, 0)) == 0);

    igraph_destroy(&pattern);
    igraph_destroy(&target);

    igraph_vector_int_destroy(&map);
    igraph_vector_int_list_destroy(&maps);


    /* Check error: pattern and target differ in directedness */
    igraph_vector_int_init(&map, 0);
    igraph_vector_int_list_init(&maps, 0);
    igraph_small(&pattern, 0, IGRAPH_DIRECTED, -1);
    CHECK_ERROR(
        igraph_subisomorphic_lad(&pattern, &target, /*domains=*/ 0,
        &iso, &map, &maps, /*induced=*/ 0),
        IGRAPH_EINVAL
    );
    igraph_vector_int_destroy(&map);
    igraph_vector_int_list_destroy(&maps);
    igraph_destroy(&pattern);
    igraph_destroy(&target);

    test_motifs();
    test_motifs_undirected();

    VERIFY_FINALLY_STACK();

    return 0;
}
