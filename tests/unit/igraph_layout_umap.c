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


int check_graph_twoclusters(const igraph_matrix_t *layout) {
    /* 4 vertices (0-3), 2 articulation points (4-5), 4 vertices (6-9) */
    igraph_real_t xm, ym, dx, dy, dist, xmin, xmax, ymin, ymax, distmax;
    int nerr = 0;

    xmin = xmax = ymin = ymax = 0;
    for (int i = 0; i < 12; i++) {
        xmin = fmin(xmin, MATRIX(*layout, i, 0));
        xmax = fmax(xmax, MATRIX(*layout, i, 0));
        ymin = fmin(ymin, MATRIX(*layout, i, 1));
        ymax = fmax(ymax, MATRIX(*layout, i, 1));
    }
    /* total span of the layout */
    distmax = fmax((xmax - xmin), (ymax - ymin));

    for (int iclu = 0; iclu < 8; iclu+= 7) {
        xm = 0;
        ym = 0;
        for (int i = iclu; i < iclu + 4; i++) {
            xm += MATRIX(*layout, i, 0);
            ym += MATRIX(*layout, i, 1);
        }
        xm /= 4;
        ym /= 4;
        for (int i = iclu; i < iclu + 4; i++) {
            dx = MATRIX(*layout, i, 0) - xm;
            dy = MATRIX(*layout, i, 1) - ym;
            dist = sqrt((dx * dx) + (dy * dy));

            if (dist > 0.2 * distmax) {
                printf("ERROR: UMAP cluster not compact!\n");
                printf("Vertex %d, dx: %f, dy: %f, dist: %f, distmax: %f, d/dmax: %e\n", i, dx, dy, dist, distmax, dist / distmax);
                nerr++;
            }
        }
    }

    if (nerr == 0) {
        printf("UMAP layout seems fine.\n");
    } else {
        igraph_matrix_print(layout);
    }
}


int main() {
    igraph_t graph, empty_graph;
    igraph_vector_t distances;
    igraph_matrix_t layout;
#ifdef UMAP_DEBUG
    igraph_real_t a, b;
#endif

    igraph_rng_seed(igraph_rng_default(), 42);
    igraph_small(&graph, 4, IGRAPH_UNDIRECTED,
            0,1, 0,2, 0,3, 1,2, 1,3, 2,3,
            3,4, 4,5, 5,6,
            6,7, 7,8, 6,8, 7,9, 6,9, 8,9, 7,10, 8,10, 9,10, 10,11, 9,11, 8,11, 7,11,
            -1);
    igraph_vector_init_real(&distances,
            igraph_ecount(&graph),
            0.1, 0.05, 0.12, 0.09, 0.1, 0.1,
            0.9, 0.9, 0.9,
            0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.08, 0.05, 0.1, 0.08, 0.12, 0.09, 0.11
            );

    igraph_matrix_init(&layout, 0, 0);

#ifdef UMAP_DEBUG
    IGRAPH_CHECK(igraph_i_umap_fit_ab(1, &a, &b));
    IGRAPH_CHECK(igraph_i_umap_fit_ab(0.1, &a, &b));
    IGRAPH_CHECK(igraph_i_umap_fit_ab(5, &a, &b));
#endif

    printf("layout of two clusters of vertices with 2 articulation points:\n");
    IGRAPH_ASSERT(igraph_layout_umap(&graph, &distances, &layout, -1, -1, -1) == IGRAPH_SUCCESS);
    check_graph_twoclusters(&layout);
#ifdef UMAP_DEBUG
    igraph_matrix_print(&layout);
#endif
    igraph_vector_destroy(&distances);

    printf("Same graph, no weights:\n");
    IGRAPH_ASSERT(igraph_layout_umap(&graph, NULL, &layout, 0.01, 500, 0.8) == IGRAPH_SUCCESS);
    check_graph_twoclusters(&layout);
#ifdef UMAP_DEBUG
    igraph_matrix_print(&layout);
#endif
    igraph_destroy(&graph);

    printf("Empty graph:\n");
	igraph_small(&empty_graph, 0, IGRAPH_UNDIRECTED, -1);
    IGRAPH_ASSERT(igraph_layout_umap(&empty_graph, NULL, &layout, 0.01, 500, -1) == IGRAPH_SUCCESS);
    igraph_matrix_print(&layout);
    igraph_destroy(&empty_graph);

    igraph_matrix_destroy(&layout);
    VERIFY_FINALLY_STACK();
    return 0;
}
