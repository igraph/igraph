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
#include "igraph_umap.h"
#include "test_utilities.inc"

int main() {
    igraph_t graph, empty_graph;
    igraph_vector_t distances;
    igraph_matrix_t layout;
    float a, b;

    igraph_rng_seed(igraph_rng_default(), 42);
    igraph_small(&graph, 4, IGRAPH_UNDIRECTED,
            0,1, 0,2, 0,3, 1,2, 1,3, 2,3,
            3,4, 4,5, 5,6,
            6,7, 7,8, 6,8, 7,9, 6,9, 8,9,
            -1);
    igraph_vector_init_real(&distances,
            igraph_ecount(&graph),
            0.1, 0.15, 0.12, 0.09, 0.1, 0.1,
            0.9, 0.9, 0.9,
            0.1, 0.1, 0.1, 0.1, 0.1, 0.1
            );

	igraph_small(&empty_graph, 0, IGRAPH_UNDIRECTED, -1);
    igraph_matrix_init(&layout, 0, 0);

#ifdef UMAP_DEBUG
    IGRAPH_CHECK(igraph_i_umap_fit_ab(1, &a, &b));
    IGRAPH_CHECK(igraph_i_umap_fit_ab(0.1, &a, &b));
    IGRAPH_CHECK(igraph_i_umap_fit_ab(5, &a, &b));
#endif
    

    printf("layout of two clusters of vertices with 2 articulation points:\n");
    IGRAPH_ASSERT(igraph_layout_umap(&graph, &distances, &layout, 0.01, 500) == IGRAPH_SUCCESS);
    igraph_matrix_print(&layout);

    printf("Same graph, no weights:\n");
    IGRAPH_ASSERT(igraph_layout_umap(&graph, NULL, &layout, 0.01, 500) == IGRAPH_SUCCESS);
    igraph_matrix_print(&layout);

    printf("Empty graph:\n");
    IGRAPH_ASSERT(igraph_layout_umap(&empty_graph, NULL, &layout, 0.01, 500) == IGRAPH_SUCCESS);
    igraph_matrix_print(&layout);


    igraph_destroy(&graph);
    igraph_matrix_destroy(&layout);
    VERIFY_FINALLY_STACK();
    return 0;
}
