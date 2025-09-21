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

int main(void) {
    igraph_t graph;
    igraph_vector_int_list_t paths, paths_edge;
    igraph_vector_int_t nrgeo;
    igraph_vector_t weights;
    igraph_int_t from, to;

    igraph_vector_int_list_init(&paths, 0);
    igraph_vector_int_list_init(&paths_edge, 0);
    igraph_vector_int_init(&nrgeo, 0);

    igraph_vector_init(&weights, 0);

    /* Note on the output:
     * get_all_shortest_paths functions sort their output based on the
     * last vertex only. Thus the ordering is not fully defined.
     *
     * This test does not currently canonicalize (i.e. sort)
     * the result before printing it.
     */

    printf("Singleton graph\n");
    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);

    from = 0; to = 0;

    igraph_get_all_shortest_paths(&graph, NULL, &paths, NULL, &nrgeo, from, igraph_vss_1(to), IGRAPH_ALL);

    printf("Vertex paths:\n");
    print_vector_int_list(&paths);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&paths) == VECTOR(nrgeo)[to]);

    printf("\nSingleton graph, weighted\n");
    igraph_vector_resize(&weights, igraph_ecount(&graph));
    igraph_vector_fill(&weights, 1);
    igraph_get_all_shortest_paths(&graph, &weights, &paths, NULL, &nrgeo, from, igraph_vss_1(to), IGRAPH_ALL);

    printf("Vertex paths:\n");
    print_vector_int_list(&paths);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&paths) == VECTOR(nrgeo)[to]);

    igraph_destroy(&graph);

    printf("\nNo paths\n");
    igraph_empty(&graph, 2, IGRAPH_UNDIRECTED);

    from = 0; to = 1;

    igraph_get_all_shortest_paths(&graph, NULL, &paths, &paths_edge, &nrgeo, from, igraph_vss_1(to), IGRAPH_ALL);

    printf("Vertex paths:\n");
    print_vector_int_list(&paths);
    printf("Edge paths:\n");
    print_vector_int_list(&paths_edge);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&paths) == VECTOR(nrgeo)[to]);

    igraph_destroy(&graph);

    /* This graph has multi-edges (which induce multiple paths of the
     * same length) as well as more paths of the same length between
     * vertices 0 and 4. */
    igraph_small(&graph, 0, IGRAPH_ADJ_UNDIRECTED,
                 0, 1, 1, 2, 2, 3, 3, 4, 0, 1, 2, 5, 4, 5, 1, 6, 6, 7, 3, 7,
                 -1);

    from = 0; to = 4;

    printf("\nUnweighted\n");
    igraph_get_all_shortest_paths(&graph, NULL, &paths, &paths_edge, &nrgeo, from, igraph_vss_1(to), IGRAPH_ALL);

    printf("Vertex paths:\n");
    print_vector_int_list(&paths);
    printf("Edge paths:\n");
    print_vector_int_list(&paths_edge);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&paths) == VECTOR(nrgeo)[to]);

    printf("\nWeighted, uniform weights\n");
    igraph_vector_resize(&weights, igraph_ecount(&graph));
    igraph_vector_fill(&weights, 1.5); /* constant weights */

    igraph_get_all_shortest_paths(&graph, &weights, &paths, &paths_edge, &nrgeo, from, igraph_vss_1(to), IGRAPH_ALL);

    printf("Vertex paths:\n");
    print_vector_int_list(&paths);
    printf("Edge paths:\n");
    print_vector_int_list(&paths_edge);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&paths) == VECTOR(nrgeo)[to]);

    printf("\nWeighted, multiple weighted shortest paths\n");
    VECTOR(weights)[1] = 3.0; /* create path with one more hop, but equal weighted length */
    VECTOR(weights)[4] = 2.0; /* break symmetry on pair of parallel edges */

    igraph_get_all_shortest_paths(&graph, &weights, &paths, &paths_edge, &nrgeo, from, igraph_vss_1(to), IGRAPH_ALL);

    printf("Vertex paths:\n");
    print_vector_int_list(&paths);
    printf("Edge paths:\n");
    print_vector_int_list(&paths_edge);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&paths) == VECTOR(nrgeo)[to]);

    igraph_vector_destroy(&weights);
    igraph_destroy(&graph);

    /* Graph is from https://github.com/igraph/rigraph/issues/314 */
    printf("\nWeighted, multiple weighted shortest paths, testing tolerances\n");
    igraph_small(&graph, 0, IGRAPH_UNDIRECTED,
                 0,3, 1,2, 2,5, 2,3, 2,4, 3,5, 5,6, 6,7, 1,8, 8,9, 4,10, 6,11, 8,12, 8,13, 4,14, 7,15,
                 -1);

    igraph_real_t weights_raw[] = { 1.9617537, 0.9060834, 2.2165446, 1.6251956,
                                    2.4473929, 0.5913490, 8.7093236, 2.8387330,
                                    6.1225042, 20.7217776, 6.8027218, 16.3147479,
                                    5.2605598, 6.6816853, 4.9482123, 1.8989790 };

    /* Choose carefully: If not using tolerances, the result would be incorrect
     * for starting vertices 5 and 6, but not for all other starting vertices. */
    from = 6;
    printf("From: %" IGRAPH_PRId ", to: all.\n", from);

    weights = igraph_vector_view(weights_raw, sizeof(weights_raw) / sizeof(igraph_real_t));
    igraph_get_all_shortest_paths(&graph, &weights, &paths, &paths_edge, &nrgeo, from, igraph_vss_all(), IGRAPH_ALL);

    printf("Vertex paths:\n");
    print_vector_int_list(&paths);
    printf("Edge paths:\n");
    print_vector_int_list(&paths_edge);

    printf("nrgeo: ");
    print_vector_int(&nrgeo);

    igraph_vector_int_list_destroy(&paths);
    igraph_vector_int_list_destroy(&paths_edge);

    igraph_vector_int_destroy(&nrgeo);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
