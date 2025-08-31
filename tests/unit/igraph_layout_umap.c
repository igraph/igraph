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


/* Helper function, computes layout span of a subset of vertices */
igraph_real_t get_layout_span(const igraph_matrix_t *layout, int i_start, int i_end) {
    igraph_real_t dx, dy, dz, xmin, xmax, ymin, ymax, zmin, zmax, spanmax;
    igraph_int_t ndim = igraph_matrix_ncol(layout);

    if (i_start >= i_end)
        return -1;

    xmin = xmax = MATRIX(*layout, i_start, 0);
    ymin = ymax = MATRIX(*layout, i_start, 1);
    if (ndim == 3) {
        zmin = zmax = MATRIX(*layout, i_start, 2);
    }
    for (int i = i_start + 1; i < i_end; i++) {
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
    dx = xmax - xmin;
    dy = ymax - ymin;
    spanmax = dx > dy ? dx : dy;
    if (ndim == 3) {
        dz = zmax - zmin;
        spanmax = dz > spanmax ? dz : spanmax;
    }
    return spanmax;
}


void check_graph_largeunion(
        const igraph_matrix_t *layout,
        const igraph_vector_int_t *subgraph_sizes) {

    /* sizes of the full subgraphs are 50, 70, 90 */
    igraph_real_t distlim;
    igraph_int_t nrow = igraph_matrix_nrow(layout);
    igraph_error_t err;
    int nerr = 0;
    igraph_int_t nvertices = 0;

    /* each cluster should occupy no more than say 50% of the layout */
    distlim = 0.2 * get_layout_span(layout, 0, nrow);
    for (int i = 0; i < igraph_vector_int_size(subgraph_sizes); i++) {
        err = get_layout_span(
                layout, nvertices, nvertices + VECTOR(*subgraph_sizes)[i]) > distlim;
        nvertices += VECTOR(*subgraph_sizes)[i];
        if (err) {
            nerr++;
            printf("layout of subgraph #%d too wide\n", i);
        }
    }
    if (nerr == 0) {
        printf("UMAP layout of large graph seems fine.\n");
    }
}


void check_graph_twoclusters_weights(
        const igraph_vector_t *weights,
        const igraph_vector_t *distances) {
    int nerr = 0;
    igraph_int_t i;
    igraph_real_t weight;
    igraph_int_t nc = igraph_vector_size(weights);
    igraph_int_t nd;

    if (distances == NULL) {
        for (i = 0; i < nc; i++) {
            if (VECTOR(*weights)[i] != 1.0) {
                printf("Connectivities with NULL distances should all be 1 or -1.\n");
                return;
            }
        }

        printf("Connectivities from NULL distances seem fine.\n");
    } else {
        nd = igraph_vector_size(distances);
        if (nd != nc) {
            nerr++;
            printf("Length of distances and weights must be equal.\n");
            return;
        }

        for (i = 0; i < nc; i++) {
            weight = VECTOR(*weights)[i];
            if (weight > 1.0) {
                printf("Weights cannot be >1.0, found %f", weight);
                return;
            } else if (weight < 0.0) {
                printf("Weights cannot be negative, found %f", weight);
                return;
            }
        }
        printf("Connectivities from distances seem fine.\n");
    }
    return;
}


void check_graph_twoclusters(const igraph_matrix_t *layout) {
    /* 4 vertices (0-3), 2 articulation points (4-5), 4 vertices (6-9) */
    igraph_real_t xm, ym, zm, dx, dy, dz, dist, distmax;
    igraph_int_t ndim = igraph_matrix_ncol(layout);

    int nerr = 0;

    distmax = get_layout_span(layout, 0, 12);

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
    igraph_int_t nrows = igraph_matrix_nrow(layout);
    igraph_int_t ncols = igraph_matrix_ncol(layout);
    igraph_int_t nerr = 0;

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


int main(void) {
    igraph_t graph, empty_graph, singleton_graph;
    igraph_vector_t distances, weights;
    igraph_matrix_t layout;
    igraph_t graph1, graph2, graph3;
    igraph_vector_ptr_t graph_ptr;
    igraph_vector_int_t subgraph_sizes;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_matrix_init(&layout, 0, 0);

    /* Check passing NULL as distances */
    igraph_ring(&graph, 10, IGRAPH_UNDIRECTED, false, true);
    igraph_layout_umap(&graph, &layout, false, NULL, 0.01, 100, true);
    igraph_destroy(&graph);

    /* Check simple graphs (empty, singleton) */
    igraph_small(&empty_graph, 0, IGRAPH_UNDIRECTED, -1);
    igraph_small(&singleton_graph, 1, IGRAPH_UNDIRECTED, -1);

    printf("Check error for negative min_dist.\n");
    CHECK_ERROR(igraph_layout_umap(&empty_graph, &layout, 0, NULL, -0.01, 500, 0), IGRAPH_EINVAL);

    printf("Check error for negative epochs.\n");
    CHECK_ERROR(igraph_layout_umap(&empty_graph, &layout, 0, NULL, 0.01, -1, 0), IGRAPH_EINVAL);

    printf("Check error for seed layout with wrong dimensions.\n");
    CHECK_ERROR(igraph_layout_umap(&empty_graph, &layout, 1, NULL, 0.01, 500, 0), IGRAPH_EINVAL);

    printf("Empty graph:\n");
    IGRAPH_ASSERT(igraph_layout_umap(&empty_graph, &layout, 0, NULL, 0.01, 500, 0) == IGRAPH_SUCCESS);
    igraph_matrix_print(&layout);

    printf("Singleton graph:\n");
    IGRAPH_ASSERT(igraph_layout_umap(&singleton_graph, &layout, 0, NULL, 0.01, 500, 0) == IGRAPH_SUCCESS);
    check_graph_singleton(&layout);

    printf("Singleton graph, use seed:\n");
    IGRAPH_ASSERT(igraph_layout_umap(&singleton_graph, &layout, 1, NULL, 0.01, 500, 0) == IGRAPH_SUCCESS);
    check_graph_singleton(&layout);

    igraph_destroy(&singleton_graph);
    igraph_destroy(&empty_graph);

    /* Check a small graph with two main groups of vertices */
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

    printf("layout of two clusters of vertices with 2 articulation points:\n");
    IGRAPH_ASSERT(igraph_layout_umap(&graph, &layout, 0, &distances, 0.01, 500, 0) == IGRAPH_SUCCESS);
    check_graph_twoclusters(&layout);

    printf("same graph, different epochs:\n");
    IGRAPH_ASSERT(igraph_layout_umap(&graph, &layout, 0, &distances, 0.0, 5000, 0) == IGRAPH_SUCCESS);
    check_graph_twoclusters(&layout);

    printf("Same graph, no distances:\n");
    IGRAPH_ASSERT(igraph_layout_umap(&graph, &layout, 0, NULL, 0.0, 500, 0) == IGRAPH_SUCCESS);
    check_graph_twoclusters(&layout);
    igraph_matrix_resize(&layout, 0, 0);

    printf("Same graph, 3D layout:\n");
    IGRAPH_ASSERT(igraph_layout_umap_3d(&graph, &layout, 0, NULL, 0.0, 500, 0) == IGRAPH_SUCCESS);
    check_graph_twoclusters(&layout);

    printf("Same graph, custom initial layout:\n");
    igraph_layout_random(&graph, &layout);
    IGRAPH_ASSERT(igraph_layout_umap(&graph, &layout, 1, NULL, 0.01, 500, 0) == IGRAPH_SUCCESS);
    check_graph_twoclusters(&layout);

    printf("Same graph, just compute weights from NULL distances:\n");
    igraph_vector_init(&weights, 0);
    igraph_layout_umap_compute_weights(&graph, NULL, &weights);
    check_graph_twoclusters_weights(&weights, NULL);

    printf("Same graph, just compute weights from distances:\n");
    igraph_layout_umap_compute_weights(&graph, &distances, &weights);
    check_graph_twoclusters_weights(&weights, &distances);

    printf("Same graph, precomputed weights:\n");
    IGRAPH_ASSERT(igraph_layout_umap(&graph, &layout, 0, &weights, 0.01, 500, 1) == IGRAPH_SUCCESS);
    check_graph_twoclusters(&layout);

    igraph_matrix_destroy(&layout);
    igraph_vector_destroy(&distances);
    igraph_vector_destroy(&weights);
    igraph_destroy(&graph);

    /* Check a larger graph with around 150 vertices, split in 3 main clusters */
    printf("Large disjoint union of 3 full graphs:\n");
    igraph_matrix_init(&layout, 0, 0);
    igraph_vector_int_init(&subgraph_sizes, 3);
    VECTOR(subgraph_sizes)[0] = 10;
    VECTOR(subgraph_sizes)[1] = 50;
    VECTOR(subgraph_sizes)[2] = 30;
    // Make 3 full graphs
    igraph_vector_ptr_init(&graph_ptr, 3);
    igraph_full(&graph1, VECTOR(subgraph_sizes)[0], 0, 0);
    VECTOR(graph_ptr)[0] = &graph1;
    igraph_full(&graph2, VECTOR(subgraph_sizes)[1], 0, 0);
    VECTOR(graph_ptr)[1] = &graph2;
    igraph_full(&graph3, VECTOR(subgraph_sizes)[2], 0, 0);
    VECTOR(graph_ptr)[2] = &graph3;
    // Get the disjoint union
    igraph_disjoint_union_many(&graph, &graph_ptr);
    // Call UMAP
    IGRAPH_ASSERT(igraph_layout_umap(
                &graph, &layout, 0, NULL, 0.0, 50, 0) == IGRAPH_SUCCESS);
    // Check the layout, it should have three balls
    check_graph_largeunion(&layout, &subgraph_sizes);

    // Destroy data structures
    igraph_matrix_destroy(&layout);
    igraph_destroy(&graph);
    igraph_destroy(&graph3);
    igraph_destroy(&graph2);
    igraph_destroy(&graph1);
    igraph_vector_int_destroy(&subgraph_sizes);
    igraph_vector_ptr_destroy(&graph_ptr);

    VERIFY_FINALLY_STACK();

    return 0;
}
