/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2021  The igraph development team <igraph@igraph.org>

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
void test_bug950(void) {
    /* Testing the case of weighted graphs with multiple alternate
     * paths to the same node with slightly different weights due to
     * floating point inaccuracies. */
    igraph_t g;
    igraph_vector_t eb;
    igraph_vector_t weights;
    igraph_integer_t from, to;
    igraph_integer_t no_of_edges, i;

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

    igraph_vector_init(&eb, 0);

    igraph_edge_betweenness(&g, &eb, IGRAPH_UNDIRECTED, &weights);
    print_vector(&eb);

    igraph_vector_destroy(&eb);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);
}


/* https://github.com/igraph/igraph/issues/1050 */
void test_bug1050(void) {
    /* compare cutoff = -1 with cutoff = 0 */
    igraph_t g;
    igraph_vector_t eb, eb2;
    igraph_vector_t weights;

    igraph_full(&g, 6, 0, 0);

    /* unweighted */
    igraph_vector_init(&eb, igraph_ecount(&g));
    igraph_vector_init(&eb2, igraph_ecount(&g));

    igraph_edge_betweenness_cutoff(&g, &eb, IGRAPH_UNDIRECTED, /* weights */ 0, /* cutoff */ -1);
    igraph_edge_betweenness_cutoff(&g, &eb2, IGRAPH_UNDIRECTED, /* weights */ 0, /* cutoff */ 0);

    /* results must differ */
    IGRAPH_ASSERT(! igraph_vector_all_e(&eb, &eb2));

    igraph_vector_destroy(&eb);
    igraph_vector_destroy(&eb2);

    /* weighted */
    igraph_vector_init(&eb, igraph_ecount(&g));
    igraph_vector_init(&eb2, igraph_ecount(&g));

    igraph_vector_init(&weights, igraph_ecount(&g));
    igraph_vector_fill(&weights, 1);
    VECTOR(weights)[0] = 2;

    igraph_edge_betweenness_cutoff(&g, &eb, IGRAPH_UNDIRECTED, &weights, /* cutoff */ -1);
    igraph_edge_betweenness_cutoff(&g, &eb2, IGRAPH_UNDIRECTED, &weights, /* cutoff */ 0);

    /* results must differ */
    IGRAPH_ASSERT(! igraph_vector_all_e(&eb, &eb2));

    igraph_vector_destroy(&eb);
    igraph_vector_destroy(&eb2);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);
}


int main(void) {
    igraph_t g;
    igraph_vector_t eb, eb2;
    igraph_vector_t weights;

    igraph_vector_init(&eb, 0);

    printf("Null graph\n");
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    igraph_edge_betweenness(&g, &eb, IGRAPH_UNDIRECTED, NULL);
    print_vector(&eb);
    igraph_destroy(&g);

    printf("\nEdgeless graph on 3 vertices\n");
    igraph_empty(&g, 3, IGRAPH_DIRECTED);
    igraph_edge_betweenness(&g, &eb, IGRAPH_DIRECTED, NULL);
    print_vector(&eb);
    igraph_destroy(&g);

    igraph_vector_destroy(&eb);

    printf("\nNo cutoff, undirected, unweighted\n");
    igraph_famous(&g, "zachary");
    igraph_vector_init(&eb, 0);
    igraph_edge_betweenness(&g, &eb, IGRAPH_UNDIRECTED, /*weights=*/ 0);
    print_vector(&eb);

    printf("\nNo cutoff, undirected, unit weighted\n");
    igraph_vector_init(&eb2, 0);
    igraph_vector_init(&weights, igraph_ecount(&g));
    igraph_vector_fill(&weights, 1.0);
    igraph_edge_betweenness(&g, &eb2, IGRAPH_UNDIRECTED, &weights);
    print_vector(&eb2);

    /* check that weighted and unweighted calculations give the same result */
    igraph_vector_scale(&eb2, -1);
    igraph_vector_add(&eb, &eb2);
    igraph_vector_abs(&eb);
    IGRAPH_ASSERT(igraph_vector_max(&eb) < 1e-13);

    igraph_vector_destroy(&weights);
    igraph_vector_destroy(&eb2);
    igraph_vector_destroy(&eb);
    igraph_destroy(&g);

    printf("\nSmall directed graph, unweighted\n");
    igraph_small(&g, 0, IGRAPH_DIRECTED,
                 1,0, 2,0, 0,3, 3,4, 4,5, 5,0, 5,6,
                 -1);
    igraph_vector_init(&eb, 0);
    igraph_edge_betweenness(&g, &eb, IGRAPH_DIRECTED, /* weights */ NULL);
    print_vector(&eb);
    igraph_vector_destroy(&eb);
    igraph_destroy(&g);

    printf("\nSmall undirected graph 1, unweighted, cutoff=2\n");
    igraph_small(&g, 0, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 1, 4, -1);
    igraph_vector_init(&eb, 0);
    igraph_edge_betweenness_cutoff(&g, &eb, IGRAPH_UNDIRECTED, /*weights=*/ 0, /*cutoff=*/2);
    print_vector(&eb);
    igraph_vector_destroy(&eb);
    igraph_destroy(&g);

    printf("\nSmall undirected graph 2, unweighted, cutoff=2\n");
    igraph_small(&g, 0, IGRAPH_UNDIRECTED,
                 0, 1, 0, 3, 1, 2, 1, 4, 2, 5, 3, 4, 3, 6, 4, 5, 4, 7, 5, 8,
                 6, 7, 7, 8, -1);
    igraph_vector_init(&eb, 0);
    igraph_edge_betweenness_cutoff(&g, &eb, IGRAPH_UNDIRECTED, /*weights=*/ 0, /*cutoff=*/2);
    print_vector(&eb);
    igraph_vector_destroy(&eb);
    igraph_destroy(&g);

    printf("\nSmall undirected graph 3, unweighted, with multiple and loop edges\n");
    igraph_small(&g, 4, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 1, 2, 1, 1, 2, 3, 3, 0, 3, 3, -1);
    igraph_vector_init(&eb, 0);
    igraph_edge_betweenness(/* graph=     */ &g,
            /* res=       */ &eb,
            /* directed = */ IGRAPH_UNDIRECTED,
            /* weights=   */ 0);
    print_vector(&eb);
    igraph_vector_destroy(&eb);
    igraph_destroy(&g);

    printf("\nTesting bug 950, tolerances\n");
    test_bug950();

    printf("\nTesting bug 1050, cutoff values\n");
    test_bug1050();

    VERIFY_FINALLY_STACK();

    return 0;
}
