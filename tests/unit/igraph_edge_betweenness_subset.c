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
#include "test_utilities.inc"


/* https://github.com/igraph/igraph/issues/950 */
void test_bug950_edge() {
    /* Testing the case of weighted graphs with multiple alternate
     * paths to the same node with slightly different weights due to
     * floating point inaccuracies. */
    igraph_t g;
    igraph_vector_t eb;
    igraph_vector_t weights;
    igraph_integer_t from, to;
    long int no_of_edges, i;

    igraph_full(&g, 6, 0, 0);
    no_of_edges = igraph_ecount(&g);

    igraph_vector_init(&weights, no_of_edges);

    for (i = 0; i < no_of_edges; i++) {
        igraph_edge(&g, i, &from, &to);
        if((from < 3 && to < 3) || (from >= 3 && to >= 3))
            VECTOR(weights)[i] = 1;
        else
            VECTOR(weights)[i] = 0.1;
    }

    printf("\nTesting bug 950 source and targets = vss_all()\n");
    printf("==========================================================\n");
    igraph_vector_init(&eb, 0);
            
    igraph_edge_betweenness_subset(/* graph=     */ &g,
        /* res=       */ &eb,
        /* eids=      */ igraph_ess_all(IGRAPH_EDGEORDER_ID),
        /* directed = */ IGRAPH_UNDIRECTED,
        /* sources = */ igraph_vss_all(),
        /* target = */ igraph_vss_all(),
        /* weights=   */ &weights);

    print_vector(&eb);
    igraph_vector_destroy(&eb);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);
}

int main() {
    igraph_t g;
    igraph_es_t es;
    igraph_vector_t eb, node_vec, source_vec, target_vec, bet, weights, edges;
    igraph_vs_t vs, vs_source, vs_target;

    /* edge betweenness test */

    {
        /* We use igraph_create() instead of igraph_small() as some MSVC versions
           will choke on an overlong argument list with "internal error C1001". */
        igraph_real_t edge_array[] = {
            0,  1,  0,  2,  0,  3,  0,  4,  0,  5,
            0,  6,  0,  7,  0,  8,  0, 10,  0, 11,
            0, 12,  0, 13,  0, 17,  0, 19,  0, 21,
            0, 31,  1,  2,  1,  3,  1,  7,  1, 13,
            1, 17,  1, 19,  1, 21,  1, 30,  2,  3,
            2,  7,  2,  8,  2,  9,  2, 13,  2, 27,
            2, 28,  2, 32,  3,  7,  3, 12,  3, 13,
            4,  6,  4, 10,  5,  6,  5, 10,  5, 16,
            6, 16,  8, 30,  8, 32,  8, 33,  9, 33,
            13, 33, 14, 32, 14, 33, 15, 32, 15, 33,
            18, 32, 18, 33, 19, 33, 20, 32, 20, 33,
            22, 32, 22, 33, 23, 25, 23, 27, 23, 29,
            23, 32, 23, 33, 24, 25, 24, 27, 24, 31,
            25, 31, 26, 29, 26, 33, 27, 33, 28, 31,
            28, 33, 29, 32, 29, 33, 30, 32, 30, 33,
            31, 32, 31, 33, 32, 33
        };
        printf("Large unweighted graph, edge betweenness\n");
        printf("==========================================================\n");
        igraph_create(&g, igraph_vector_view(&edges, edge_array, sizeof(edge_array) / sizeof(igraph_real_t)), 0, IGRAPH_UNDIRECTED);
        igraph_vector_init_seq(&source_vec, 0, 32);
        igraph_vs_vector(&vs_source, &source_vec);
        igraph_vector_init_seq(&target_vec, 1, 33);
        igraph_vs_vector(&vs_target, &target_vec);
        igraph_vector_init(&eb, 0);

        igraph_edge_betweenness_subset (/* graph=     */ &g,
        /* res=       */ &eb,
        /* eids=      */ igraph_ess_all(IGRAPH_EDGEORDER_ID),
        /* directed = */ IGRAPH_UNDIRECTED,
        /* sources = */ vs_source,
        /* target = */ vs_target,
        /* weights=   */ NULL);

        print_vector(&eb);
        igraph_vector_destroy(&eb);
        igraph_destroy(&g);
        igraph_vs_destroy(&vs_source);
        igraph_vector_destroy(&source_vec);
        igraph_vs_destroy(&vs_target);
        igraph_vector_destroy(&target_vec);
    }
        
    printf("\nSmall unweighted graph, edge betweenness\n");
    printf("==========================================================\n");
    igraph_small(&g, 0, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 1, 4, -1);

    igraph_vector_init(&eb, 0);
    igraph_vector_init_seq(&node_vec, 0, 3);
    igraph_es_vector(&es, &node_vec);
    igraph_vector_init_seq(&target_vec, 1, 4);
    igraph_vs_vector(&vs_target, &target_vec);
    igraph_edge_betweenness_subset (/* graph=     */ &g,
        /* res=       */ &eb,
        /* eids=      */ es,
        /* directed = */ IGRAPH_UNDIRECTED,
        /* sources = */ igraph_vss_all(),
        /* target = */ vs_target,
        /* weights=   */ NULL);

    print_vector(&eb);
    igraph_vector_destroy(&eb);
    igraph_es_destroy(&es);
    igraph_vector_destroy(&node_vec);
    igraph_vs_destroy(&vs_target);
    igraph_vector_destroy(&target_vec);
    igraph_destroy(&g);

    printf("\nSmall unweighted graph, edge betweenness with subset of sources\n");
    printf("==========================================================\n");
    igraph_small(&g, 0, IGRAPH_UNDIRECTED,
                 0, 1, 0, 3, 1, 2, 1, 4, 2, 5, 3, 4, 3, 6, 4, 5, 4, 7, 5, 8,
                 6, 7, 7, 8, -1);
    igraph_vector_init(&eb, 0);
    igraph_vector_init_seq(&source_vec, 1, 8);
    igraph_vs_vector(&vs_source, &source_vec);

    igraph_edge_betweenness_subset (/* graph=     */ &g,
        /* res=       */ &eb,
        /* eids=      */ igraph_ess_all(IGRAPH_EDGEORDER_ID),
        /* directed = */ IGRAPH_UNDIRECTED,
        /* sources = */ vs_source,
        /* target = */ igraph_vss_all(),
        /* weights=   */ NULL);
    print_vector(&eb);
    igraph_vector_destroy(&eb);
    igraph_vs_destroy(&vs_source);
    igraph_vector_destroy(&source_vec);
    igraph_destroy(&g);

    test_bug950_edge();

    printf("\nEmpty graph\n");
    printf("==========================================================\n");
    igraph_empty(&g, 2, IGRAPH_UNDIRECTED);
    
    igraph_vector_init(&bet, 0);
    igraph_edge_betweenness_subset (/* graph=     */ &g,
        /* res=       */ &bet,
        /* eids=      */ igraph_ess_all(IGRAPH_EDGEORDER_ID),
        /* directed = */ IGRAPH_UNDIRECTED,
        /* sources = */ igraph_vss_all(),
        /* target = */ igraph_vss_all(),
        /* weights=   */ NULL);
    print_vector(&bet);
    igraph_vector_destroy(&bet);

    igraph_destroy(&g);

    printf("\n37x37 grid graph\n");
    printf("==========================================================\n");

    {
        igraph_vector_t dims;

        igraph_vector_init(&dims, 2);
        VECTOR(dims)[0] = 37;
        VECTOR(dims)[1] = 37;

        igraph_lattice(&g, &dims, 1, IGRAPH_UNDIRECTED, 0, 0);

        igraph_vector_init(&bet, 0);
        igraph_vector_init_seq(&target_vec, 0, (int) igraph_vcount(&g) - 1);
        igraph_vector_remove(&target_vec, (long int)0);
        igraph_vs_vector(&vs_target, &target_vec);
        igraph_vector_init_seq(&source_vec, 0, (int) igraph_vcount(&g) - 1);
        igraph_vector_remove(&source_vec, (long int) 0);
        igraph_vs_vector(&vs_source, &source_vec);
        
        igraph_edge_betweenness_subset (/* graph=     */ &g,
        /* res=       */ &bet,
        /* eids=      */ igraph_ess_all(IGRAPH_EDGEORDER_ID),
        /* directed = */ IGRAPH_UNDIRECTED,
        /* sources = */ vs_source,
        /* target = */ vs_target,
        /* weights=   */ NULL);
        printf("Max edge betweenness: %f\n", igraph_vector_max(&bet));

        igraph_vector_destroy(&bet);
        igraph_destroy(&g);
        igraph_vector_destroy(&dims);
        igraph_vs_destroy(&vs_target);
        igraph_vector_destroy(&target_vec);
        igraph_vs_destroy(&vs_source);
        igraph_vector_destroy(&source_vec);
    }

    VERIFY_FINALLY_STACK();

    return 0;
}
