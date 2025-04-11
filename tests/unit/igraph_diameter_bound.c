/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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

int main(void) {
    igraph_t g;
    igraph_vector_t weights;
    igraph_real_t result;
    igraph_real_t reference;

    // recreate example graph from paper. A=0, B=1, etc. Note that O is missing.
    // result should be 7
    printf("Sample graph with diameter 7\n");
    igraph_small(&g, 19, IGRAPH_UNDIRECTED, 
        0,2, 1,2, 2,3, 2,4, 2,5, 3,5, 4,5, 4,6, 5,6, 5,7, 5,9, 6,8, 6,9, 7,10,
        9,11, 11,12, 11,13, 11,14, 13,14, 14,15, 14,16, 15,17, 16,18, -1
    );
    igraph_diameter_bound(&g, NULL, &result, IGRAPH_UNDIRECTED, 0);
    igraph_diameter(&g, &reference, NULL, NULL, NULL, NULL, IGRAPH_UNDIRECTED, 1);
    IGRAPH_ASSERT(result == reference && result == 7);
    igraph_destroy(&g);

    // ring graph - diameter should be half of the node count
    printf("\nRing graph with 24 nodes\n");
    igraph_ring(&g, 24, IGRAPH_UNDIRECTED, false, true);
    igraph_diameter_bound(&g, NULL, &result, IGRAPH_UNDIRECTED, 0);
    IGRAPH_ASSERT(result == 12);
    igraph_destroy(&g);

    // weighted ring graph
    igraph_vector_init_real(&weights, 9, 1.0, 2.675, 3.0, 4.0, 5.5, 1.0, 1.0, 1.0, 1.0);
    igraph_ring(&g, 10, IGRAPH_UNDIRECTED, 0, 0);
    igraph_diameter_bound(&g, &weights, &result, false, true);
    igraph_diameter_dijkstra(&g, &weights, &reference, NULL, NULL, NULL, NULL, false, true);
    IGRAPH_ASSERT(igraph_almost_equals(result, reference, 1e-10));
    printf("weighted results: %f\n%f\n\n\n", result, reference);
    igraph_destroy(&g);

    // random barabasi - compare to default diameter. Do it a few times
    printf("\nBarabasi random graph\n");
    for (igraph_integer_t i = 0; i < 10; i++) {
        igraph_rng_seed(igraph_rng_default(), 1234*i);
        igraph_barabasi_game(&g, 100, 1, 4, 0, 0, 1, IGRAPH_UNDIRECTED, IGRAPH_BARABASI_BAG, NULL);
        igraph_diameter_bound(&g, NULL, &result, IGRAPH_UNDIRECTED, 0);
        igraph_diameter(&g, &reference, NULL, NULL, NULL, NULL, IGRAPH_UNDIRECTED, 1);
        IGRAPH_ASSERT(result == reference);  // no reference value to compare to
        igraph_destroy(&g);
    }

    // random watts strogatz - compare to default diameter. Do it a few times
    printf("\nWatts Strogatz random graph\n");
    for (igraph_integer_t i = 0; i < 10; i++) {
        igraph_rng_seed(igraph_rng_default(), 1234*i);
        igraph_watts_strogatz_game(&g, 1, 100, 5, 0.2, false, false);
        igraph_diameter_bound(&g, NULL, &result, IGRAPH_UNDIRECTED, 0);
        igraph_diameter(&g, &reference, NULL, NULL, NULL, NULL, IGRAPH_UNDIRECTED, 1);
        IGRAPH_ASSERT(result == reference);  // no reference value to compare to
        igraph_destroy(&g);
    }

    //disconnected graph
    igraph_vector_int_t edge_vec;
    igraph_integer_t vec[] = {2,8};
    igraph_es_t edge_sele;
    printf("\nDisconnected ring graph, unconn=false\n");
    igraph_ring(&g, 10, IGRAPH_UNDIRECTED, 0, 0);
    igraph_vector_int_view(&edge_vec, vec, sizeof(vec) / sizeof(vec[0]));
    igraph_es_vector(&edge_sele, &edge_vec);
    igraph_delete_edges(&g, edge_sele);
    igraph_diameter_bound(&g, NULL, &result, IGRAPH_UNDIRECTED, false);
    igraph_diameter(&g, &reference, NULL, NULL, NULL, NULL, IGRAPH_DIRECTED, 0);
    IGRAPH_ASSERT(result == reference);
    printf("\nDisconnected ring graph, unconn=true\n");
    igraph_diameter_bound(&g, NULL, &result, IGRAPH_UNDIRECTED, true);
    igraph_diameter(&g, &reference, NULL, NULL, NULL, NULL, IGRAPH_DIRECTED, 1);
    IGRAPH_ASSERT(result == reference);
    igraph_es_destroy(&edge_sele);
    igraph_destroy(&g);

    // graph with zero nodes - result should be NaN
    printf("\nGraph with no vertices\n");
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    igraph_diameter_bound(&g, NULL, &result, IGRAPH_UNDIRECTED, 0);
    IGRAPH_ASSERT(result != result);  // NaN test
    igraph_destroy(&g);

    // graph with one node - result should be 0
    printf("\nGraph with one vertex\n");
    igraph_empty(&g, 1, IGRAPH_UNDIRECTED);
    igraph_diameter_bound(&g, NULL, &result, IGRAPH_UNDIRECTED, 0);
    igraph_diameter(&g, &reference, NULL, NULL, NULL, NULL, IGRAPH_UNDIRECTED, 0);
    IGRAPH_ASSERT(result == reference && result == 0);
    igraph_destroy(&g);

    igraph_vector_destroy(&weights);
    VERIFY_FINALLY_STACK();
    return 0;
}
