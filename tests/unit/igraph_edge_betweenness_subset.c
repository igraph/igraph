/*
   igraph library.
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


/* https://github.com/igraph/igraph/issues/950 */
void test_bug950_edge(void) {
    /* Testing the case of weighted graphs with multiple alternate
     * paths to the same node with slightly different weights due to
     * floating point inaccuracies. */
    igraph_t g;
    igraph_vector_t eb;
    igraph_vector_t weights;
    igraph_int_t from, to;
    igraph_int_t no_of_edges, i;

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
            /* weights=   */ &weights,
            /* res=       */ &eb,
            /* sources = */ igraph_vss_all(),
            /* target = */ igraph_vss_all(),
            /* eids=      */ igraph_ess_all(IGRAPH_EDGEORDER_ID),
            /* directed = */ IGRAPH_UNDIRECTED, false);

    print_vector(&eb);
    igraph_vector_destroy(&eb);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);
}

int main(void) {
    igraph_t g;
    igraph_es_t es;
    igraph_vector_t eb, bet, bet2, weights;
    igraph_vector_int_t node_vec, source_vec, target_vec;
    igraph_vs_t vs_source, vs_target;
    igraph_int_t i, n;

    /* edge betweenness test */

    printf("Tree\n");
    printf("==========================================================\n");
    igraph_kary_tree(&g, 11111, 10, IGRAPH_TREE_UNDIRECTED);

    /* We are including the rightmost 200 vertices from the lowermost layer
     * (layer 5) of the tree. These have 20 parents in layer 4, 2 grandparents
     * in layer 3, and a single grand-grandparent in layer 2. None of the
     * shortest paths we consider should therefore pass through the root, the
     * first 9 vertices of layer 2, the first 98 vertices of layer 3 or the
     * first 980 vertices of layer 4; the edge betweenness of the "upward"
     * pointing edges from these vertices should all be zeros.
     *
     * (Note that due to how igraph constructs tree graphs, edge i is always the
     * edge that leads from vertex i+1 towards the root).
     *
     * Also, the edge betweenness of the edges leading to children of the common
     * grand-grandparent in layer 2 is easy to calculate as any shortest path
     * going between grand-grandchildren reachable via its left child and via
     * its right child should pass through them. This gives us a betweenness of
     * 100 * 100 = 10000 for both of these edges.
     *
     * Similar calculations reveal that the edge betwennesses in subsequent
     * layers are 1900 (10 * 190) and 199 (1 * 199) if the edge being considered
     * leads to a descendant of the common grand-grandparent in layer 2, and
     * zero otherwise. */

    igraph_vs_range(&vs_source, 10911, 11111);
    igraph_vs_range(&vs_target, 10911, 11111);
    igraph_vector_init(&bet, 0);

    igraph_edge_betweenness_subset(
            /* graph=     */ &g,
            /* weights=   */ NULL,
            /* res=       */ &bet,
            /* sources =  */ vs_source,
            /* target =   */ vs_target,
            /* eids=      */ igraph_ess_all(IGRAPH_EDGEORDER_ID),
            /* directed = */ IGRAPH_UNDIRECTED, false);

    printf("Max edge betweenness: %f\n", igraph_vector_max(&bet));

    n = igraph_ecount(&g);
    for (i = 0; i < n; i++) {
        igraph_int_t expected;
        igraph_int_t vid = i + 1;

        if (vid >= 10911) {
            /* edge leading to layer 5, in the subset. There are 199 shortest
             * paths that pass through this edge, one path to any _other_
             * node in the selected subset of 200 nodes */
            expected = 199;
        } else if (vid >= 1111) {
            /* edge leading to layer 5, not in the subset */
            expected = 0;
        } else if (vid >= 1091) {
            /* edge leading to layer 4, rightmost 20 nodes */
            expected = 1900;
        } else if (vid >= 111) {
            /* edge leading to layer 4, remaining nodes */
            expected = 0;
        } else if (vid >= 109) {
            /* edge leading to layer 3, rightmost 2 nodes */
            expected = 10000;
        } else if (vid >= 11) {
            /* edge leading to layer 3, remaining nodes */
            expected = 0;
        } else if (vid == 10) {
            /* edge leading to layer 2, rightmost node */
            expected = 0;
        } else {
            expected = 0;
        }

        if (VECTOR(bet)[i] != expected) {
            printf(
                "Invalid betweenness for edge %" IGRAPH_PRId " (from vertex %" IGRAPH_PRId " towards the "
                "root), expected %" IGRAPH_PRId ", got %" IGRAPH_PRId "\n",
                i, vid, expected, (igraph_int_t) VECTOR(bet)[i]
            );
            break;
        }
    }

    igraph_vector_init(&bet2, 0);
    igraph_vector_init(&weights, igraph_ecount(&g));
    igraph_vector_fill(&weights, 1.0);

    igraph_edge_betweenness_subset(
            /* graph=     */ &g,
            /* weights=   */ &weights,
            /* res=       */ &bet2,
            /* sources = */ vs_source,
            /* target = */ vs_target,
            /* eids=      */ igraph_ess_all(IGRAPH_EDGEORDER_ID),
            /* directed = */ IGRAPH_UNDIRECTED, false);

    IGRAPH_ASSERT(igraph_vector_all_e(&bet, &bet2));

    igraph_vector_destroy(&weights);
    igraph_vs_destroy(&vs_source);
    igraph_vs_destroy(&vs_target);
    igraph_vector_destroy(&bet);
    igraph_vector_destroy(&bet2);
    igraph_destroy(&g);

    printf("\nZachary karate club, unweighted graph, edge betweenness\n");
    printf("==========================================================\n");
    igraph_famous(&g, "zachary");
    igraph_vector_int_init_range(&source_vec, 0, 33);
    igraph_vs_vector(&vs_source, &source_vec);
    igraph_vector_int_init_range(&target_vec, 1, 34);
    igraph_vs_vector(&vs_target, &target_vec);
    igraph_vector_init(&eb, 0);

    igraph_edge_betweenness_subset(/* graph=     */ &g,
            /* weights=   */ NULL,
            /* res=       */ &eb,
            /* sources = */ vs_source,
            /* target = */ vs_target,
            /* eids=      */ igraph_ess_all(IGRAPH_EDGEORDER_ID),
            /* directed = */ IGRAPH_UNDIRECTED, false);

    print_vector(&eb);
    igraph_vector_destroy(&eb);
    igraph_destroy(&g);
    igraph_vs_destroy(&vs_source);
    igraph_vector_int_destroy(&source_vec);
    igraph_vs_destroy(&vs_target);
    igraph_vector_int_destroy(&target_vec);

    printf("\nSmall unweighted graph, edge betweenness\n");
    printf("==========================================================\n");
    igraph_small(&g, 0, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 1, 4, -1);

    igraph_vector_init(&eb, 0);
    igraph_vector_int_init_range(&node_vec, 0, 4);
    igraph_es_vector(&es, &node_vec);
    igraph_vector_int_init_range(&target_vec, 1, 5);
    igraph_vs_vector(&vs_target, &target_vec);
    igraph_edge_betweenness_subset(/* graph=     */ &g,
            /* weights=   */ NULL,
            /* res=       */ &eb,
            /* sources = */ igraph_vss_all(),
            /* target = */ vs_target,
            /* eids=      */ es,
            /* directed = */ IGRAPH_UNDIRECTED, false);

    print_vector(&eb);
    igraph_vector_destroy(&eb);
    igraph_es_destroy(&es);
    igraph_vector_int_destroy(&node_vec);
    igraph_vs_destroy(&vs_target);
    igraph_vector_int_destroy(&target_vec);
    igraph_destroy(&g);

    printf("\nSmall unweighted graph, edge betweenness with subset of sources\n");
    printf("==========================================================\n");
    igraph_small(&g, 0, IGRAPH_UNDIRECTED,
                 0, 1, 0, 3, 1, 2, 1, 4, 2, 5, 3, 4, 3, 6, 4, 5, 4, 7, 5, 8,
                 6, 7, 7, 8, -1);
    igraph_vector_init(&eb, 0);
    igraph_vector_int_init_range(&source_vec, 1, 9);
    igraph_vs_vector(&vs_source, &source_vec);

    igraph_edge_betweenness_subset(/* graph=     */ &g,
            /* weights=   */ NULL,
            /* res=       */ &eb,
            /* sources = */ vs_source,
            /* target = */ igraph_vss_all(),
            /* eids=      */ igraph_ess_all(IGRAPH_EDGEORDER_ID),
            /* directed = */ IGRAPH_UNDIRECTED, false);
    print_vector(&eb);
    igraph_vector_destroy(&eb);
    igraph_vs_destroy(&vs_source);
    igraph_vector_int_destroy(&source_vec);
    igraph_destroy(&g);

    test_bug950_edge();

    printf("\nEmpty graph\n");
    printf("==========================================================\n");
    igraph_empty(&g, 2, IGRAPH_UNDIRECTED);

    igraph_vector_init(&bet, 0);
    igraph_edge_betweenness_subset(/* graph=     */ &g,
            /* weights=   */ NULL,
            /* res=       */ &bet,
            /* sources = */ igraph_vss_all(),
            /* target = */ igraph_vss_all(),
            /* eids=      */ igraph_ess_all(IGRAPH_EDGEORDER_ID),
            /* directed = */ IGRAPH_UNDIRECTED, false);
    print_vector(&bet);
    igraph_vector_destroy(&bet);

    igraph_destroy(&g);

    printf("\n37x37 grid graph\n");
    printf("==========================================================\n");

    {
        igraph_vector_int_t dims;

        igraph_vector_int_init(&dims, 2);
        VECTOR(dims)[0] = 37;
        VECTOR(dims)[1] = 37;

        igraph_square_lattice(&g, &dims, 1, IGRAPH_UNDIRECTED, /* mutual */ 0, /* periodic */ 0);

        igraph_vector_init(&bet, 0);
        igraph_vector_int_init_range(&target_vec, 0, igraph_vcount(&g));
        igraph_vector_int_remove(&target_vec, 0);
        igraph_vs_vector(&vs_target, &target_vec);
        igraph_vector_int_init_range(&source_vec, 0, igraph_vcount(&g));
        igraph_vector_int_remove(&source_vec, 0);
        igraph_vs_vector(&vs_source, &source_vec);

        igraph_edge_betweenness_subset(/* graph=     */ &g,
                /* weights=   */ NULL,
                /* res=       */ &bet,
                /* sources = */ vs_source,
                /* target = */ vs_target,
                /* eids=      */ igraph_ess_all(IGRAPH_EDGEORDER_ID),
                /* directed = */ IGRAPH_UNDIRECTED, false);
        printf("Max edge betweenness: %f\n", igraph_vector_max(&bet));

        igraph_vector_destroy(&bet);
        igraph_destroy(&g);
        igraph_vector_int_destroy(&dims);
        igraph_vs_destroy(&vs_target);
        igraph_vector_int_destroy(&target_vec);
        igraph_vs_destroy(&vs_source);
        igraph_vector_int_destroy(&source_vec);
    }

    VERIFY_FINALLY_STACK();

    return 0;
}
