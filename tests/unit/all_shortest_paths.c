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

int main() {
    igraph_t graph;
    igraph_vector_ptr_t paths;
    igraph_vector_t nrgeo, weights;
    igraph_integer_t from, to;
    long int i;

    igraph_vector_ptr_init(&paths, 0);
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&paths, igraph_vector_destroy);
    igraph_vector_init(&nrgeo, 0);

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

    igraph_get_all_shortest_paths(&graph, &paths, &nrgeo, from, igraph_vss_1(to), IGRAPH_ALL);

    for (i=0; i < igraph_vector_ptr_size(&paths); ++i) {
        print_vector(VECTOR(paths)[i]);
    }
    IGRAPH_ASSERT(igraph_vector_ptr_size(&paths) == VECTOR(nrgeo)[to]);

    igraph_vector_ptr_free_all(&paths);

    printf("\nSingleton graph, weighted\n");
    igraph_vector_resize(&weights, igraph_ecount(&graph));
    igraph_vector_fill(&weights, 1);
    igraph_get_all_shortest_paths_dijkstra(&graph, &paths, &nrgeo, from, igraph_vss_1(to), &weights, IGRAPH_ALL);

    for (i=0; i < igraph_vector_ptr_size(&paths); ++i) {
        print_vector(VECTOR(paths)[i]);
    }
    IGRAPH_ASSERT(igraph_vector_ptr_size(&paths) == VECTOR(nrgeo)[to]);

    igraph_vector_ptr_free_all(&paths);

    igraph_destroy(&graph);

    printf("\nNo paths\n");
    igraph_empty(&graph, 2, IGRAPH_UNDIRECTED);

    from = 0; to = 1;

    igraph_get_all_shortest_paths(&graph, &paths, &nrgeo, from, igraph_vss_1(to), IGRAPH_ALL);

    for (i=0; i < igraph_vector_ptr_size(&paths); ++i) {
        print_vector(VECTOR(paths)[i]);
    }
    IGRAPH_ASSERT(igraph_vector_ptr_size(&paths) == VECTOR(nrgeo)[to]);

    igraph_vector_ptr_free_all(&paths);

    igraph_destroy(&graph);

    /* This graph has multi-edges (which induce multiple paths of the
     * same length) as well as more paths of the same length between
     * vertices 0 and 4. */
    igraph_small(&graph, 0, IGRAPH_ADJ_UNDIRECTED,
                 0, 1, 1, 2, 2, 3, 3, 4, 0, 1, 2, 5, 4, 5, 1, 6, 6, 7, 3, 7,
                 -1);

    from = 0; to = 4;

    printf("\nUnweighted\n");
    igraph_get_all_shortest_paths(&graph, &paths, &nrgeo, from, igraph_vss_1(to), IGRAPH_ALL);

    for (i=0; i < igraph_vector_ptr_size(&paths); ++i) {
        print_vector(VECTOR(paths)[i]);
    }
    IGRAPH_ASSERT(igraph_vector_ptr_size(&paths) == VECTOR(nrgeo)[to]);

    igraph_vector_ptr_free_all(&paths);

    printf("\nWeighted, uniform weights\n");
    igraph_vector_resize(&weights, igraph_ecount(&graph));
    igraph_vector_fill(&weights, 1.5); /* constant weights */

    igraph_get_all_shortest_paths_dijkstra(&graph, &paths, &nrgeo, from, igraph_vss_1(to), &weights, IGRAPH_ALL);

    for (i=0; i < igraph_vector_ptr_size(&paths); ++i) {
        print_vector(VECTOR(paths)[i]);
    }
    IGRAPH_ASSERT(igraph_vector_ptr_size(&paths) == VECTOR(nrgeo)[to]);

    igraph_vector_ptr_free_all(&paths);

    printf("\nWeighted, multiple weighted shortest paths\n");
    VECTOR(weights)[1] = 3.0; /* create path with one more hop, but equal weighted length */
    VECTOR(weights)[4] = 2.0; /* break symmetry on pair of parallel edges */

    igraph_get_all_shortest_paths_dijkstra(&graph, &paths, &nrgeo, from, igraph_vss_1(to), &weights, IGRAPH_ALL);

    for (i=0; i < igraph_vector_ptr_size(&paths); ++i) {
        print_vector(VECTOR(paths)[i]);
    }
    IGRAPH_ASSERT(igraph_vector_ptr_size(&paths) == VECTOR(nrgeo)[to]);

    igraph_vector_ptr_destroy_all(&paths);

    igraph_vector_destroy(&weights);
    igraph_vector_destroy(&nrgeo);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
