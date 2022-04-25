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


void check_graph_twoclusters(const igraph_matrix_t *layout) {
    /* 4 vertices (0-3), 2 articulation points (4-5), 4 vertices (6-9) */
    igraph_real_t xm, ym, zm, dx, dy, dz, dist, xmin, xmax, ymin, ymax, zmin, zmax, distmax;
    igraph_integer_t ndim = igraph_matrix_ncol(layout);
    int nerr = 0;

    xmin = xmax = ymin = ymax = zmin = zmax = 0;
    for (int i = 0; i < 12; i++) {
        xmin = fmin(xmin, MATRIX(*layout, i, 0));
        xmax = fmax(xmax, MATRIX(*layout, i, 0));
        ymin = fmin(ymin, MATRIX(*layout, i, 1));
        ymax = fmax(ymax, MATRIX(*layout, i, 1));
        if (ndim == 3) {
            zmin = fmin(zmin, MATRIX(*layout, i, 2));
            zmax = fmax(zmax, MATRIX(*layout, i, 2));
        }
    }
    /* total span of the layout */
    distmax = fmax((xmax - xmin), (ymax - ymin));
    if (ndim == 3) {
        distmax = fmax((zmax - zmin), distmax);
    }

    for (int iclu = 0; iclu < 8; iclu+= 7) {
        xm = 0;
        ym = 0;
        zm = 0;
        for (int i = iclu; i < iclu + 4; i++) {
            xm += MATRIX(*layout, i, 0);
            ym += MATRIX(*layout, i, 1);
            if (ndim == 3) {
                zm += MATRIX(*layout, i, 2);
            }
        }
        xm /= 4;
        ym /= 4;
        zm /= 4;
        for (int i = iclu; i < iclu + 4; i++) {
            dx = MATRIX(*layout, i, 0) - xm;
            dy = MATRIX(*layout, i, 1) - ym;
            if (ndim == 3) {
                dz = MATRIX(*layout, i, 2) - zm;
            }
            dist = (dx * dx) + (dy * dy);
            if (ndim == 3) {
                dist += (dz * dz);
            }
            dist = sqrt(dist);

            if (dist > 0.2 * distmax) {
                printf("ERROR: UMAP cluster not compact!\n");
                printf("Cluster %d of 2, vertex %d, dx: %g, dy: %g, dist: %g, distmax: %g, d/dmax: %e\n", 1 + iclu / 7, i, dx, dy, dist, distmax, dist / distmax);
                nerr++;
            }
        }
    }

    if (nerr == 0) {
        printf("UMAP layout seems fine.\n");
    } else {
        printf("UMAP layout inconsistent:\nx\ty\n");
        igraph_matrix_print(layout);
    }
}


void check_graph_singleton(const igraph_matrix_t *layout) {
    igraph_integer_t nrows = igraph_matrix_nrow(layout);
    igraph_integer_t ncols = igraph_matrix_ncol(layout);
    igraph_integer_t nerr = 0;

    if (nrows != 1) {
        printf("Singleton graph layout has %d rows instead of 1.\n", (int)nrows);
        nerr++;
    }
    if (ncols != 2) {
        printf("Singleton graph layout has %d cols instead of 2.\n", (int)ncols);
        nerr++;
    }
    if ((fabs(MATRIX(*layout, 0, 0)) > 0.001) || (fabs(MATRIX(*layout, 0, 1)) > 0.001)) {
        printf("Singleton graph layout is not (0,0): (%g,%g)\n", MATRIX(*layout, 0, 0), MATRIX(*layout, 0, 1));
        nerr++;
    }
    if (nerr == 0) {
        printf("UMAP layout seems fine.\n");
    }
}


int main() {
    igraph_t graph, empty_graph, singleton_graph;
    igraph_vector_t distances;
    igraph_matrix_t layout;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_matrix_init(&layout, 0, 0);

    igraph_small(&empty_graph, 0, IGRAPH_UNDIRECTED, -1);
    igraph_small(&singleton_graph, 1, IGRAPH_UNDIRECTED, -1);
    igraph_small(&graph, 12, IGRAPH_UNDIRECTED,
            0,1, 0,2, 0,3, 1,2, 1,3, 2,3,
            3,4, 4,5, 5,6,
            6,7, 7,8, 6,8, 7,9, 6,9, 8,9, 7,10, 8,10, 9,10, 10,11, 9,11, 8,11, 7,11,
            -1);
    igraph_vector_init_real(&distances,
            igraph_ecount(&graph),
            0.1, 0.09, 0.12, 0.09, 0.1, 0.1,
            0.9, 0.9, 0.9,
            0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.08, 0.05, 0.1, 0.08, 0.12, 0.09, 0.11
            );

    printf("Check error for negative min_dist.\n");
    CHECK_ERROR(igraph_layout_umap(&empty_graph, &layout, 0, NULL, -0.01, 500, 1), IGRAPH_EINVAL);

    printf("Check error for negative epochs.\n");
    CHECK_ERROR(igraph_layout_umap(&empty_graph, &layout, 0, NULL, 0.01, -1, 1), IGRAPH_EINVAL);

    printf("Check error for negative sampling probability.\n");
    CHECK_ERROR(igraph_layout_umap(&empty_graph, &layout, 0, NULL, 0.01, 500, -1), IGRAPH_EINVAL);

    printf("Check error for sampling probability above one.\n");
    CHECK_ERROR(igraph_layout_umap(&empty_graph, &layout, 0, NULL, 0.01, 500, 1.1), IGRAPH_EINVAL);

    printf("Check error for seed layout with wrong dimensions.\n");
    CHECK_ERROR(igraph_layout_umap(&empty_graph, &layout, 1, NULL, 0.01, 500, 1), IGRAPH_EINVAL);

    printf("Empty graph:\n");
    IGRAPH_ASSERT(igraph_layout_umap(&empty_graph, &layout, 0, NULL, 0.01, 500, 1) == IGRAPH_SUCCESS);
    igraph_matrix_print(&layout);
    igraph_destroy(&empty_graph);

    printf("Singleton graph:\n");
    IGRAPH_ASSERT(igraph_layout_umap(&singleton_graph, &layout, 0, NULL, 0.01, 500, 0.2) == IGRAPH_SUCCESS);
    check_graph_singleton(&layout);
    igraph_destroy(&singleton_graph);

    printf("layout of two clusters of vertices with 2 articulation points:\n");
    IGRAPH_ASSERT(igraph_layout_umap(&graph, &layout, 0, &distances, 0.01, 500, 0.3) == IGRAPH_SUCCESS);
    check_graph_twoclusters(&layout);

    printf("same graph, different negative sampling probability:\n");
    IGRAPH_ASSERT(igraph_layout_umap(&graph, &layout, 0, &distances, 0.01, 500, 0.8) == IGRAPH_SUCCESS);
    check_graph_twoclusters(&layout);

    printf("same graph, different epochs:\n");
    IGRAPH_ASSERT(igraph_layout_umap(&graph, &layout, 0, &distances, 0.01, 5000, 0.8) == IGRAPH_SUCCESS);
    check_graph_twoclusters(&layout);
    igraph_vector_destroy(&distances);

    printf("Same graph, no distances:\n");
    IGRAPH_ASSERT(igraph_layout_umap(&graph, &layout, 0, NULL, 0.01, 500, 0.8) == IGRAPH_SUCCESS);
    check_graph_twoclusters(&layout);
    igraph_matrix_resize(&layout, 0, 0);

    printf("Same graph, 3D layout:\n");
    IGRAPH_ASSERT(igraph_layout_umap_3d(&graph, &layout, 0, NULL, 0.01, 500, 0.8) == IGRAPH_SUCCESS);
    check_graph_twoclusters(&layout);

    printf("Same graph, custom initial layout:\n");
    igraph_layout_random(&graph, &layout);
    IGRAPH_ASSERT(igraph_layout_umap(&graph, &layout, 1, NULL, 0.01, 500, 0.8) == IGRAPH_SUCCESS);
    check_graph_twoclusters(&layout);

    igraph_matrix_destroy(&layout);
    igraph_destroy(&graph);
    VERIFY_FINALLY_STACK();

    return 0;
}
