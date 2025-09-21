/*
   igraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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
#include <math.h>

#include "test_utilities.h"

void check_radius(const igraph_t *graph, const igraph_vector_int_t *center, igraph_neimode_t mode) {
    igraph_vector_t ecc;
    igraph_real_t radius;
    igraph_int_t n = igraph_vector_int_size(center);

    igraph_radius(graph, NULL, &radius, mode);
    printf("Radius: %g\n", radius);

    if (n == 0) {
        /* Null graph has radius NaN */
        IGRAPH_ASSERT(isnan(radius));
    } else {
        igraph_vector_init(&ecc, 0);
        igraph_eccentricity(graph, NULL, &ecc, igraph_vss_vector(center), mode);
        for (igraph_int_t i=0; i < n; i++) {
            IGRAPH_ASSERT(VECTOR(ecc)[i] == radius);
        }
        igraph_vector_destroy(&ecc);
    }
}

void check_radius_dijkstra(const igraph_t *graph, const igraph_vector_t *weights,
                           const igraph_vector_int_t *center, igraph_neimode_t mode) {
    igraph_vector_t ecc;
    igraph_real_t radius;
    igraph_int_t n = igraph_vector_int_size(center);
    const igraph_real_t eps = IGRAPH_SHORTEST_PATH_EPSILON;

    igraph_radius(graph, weights, &radius, mode);
    printf("Radius: %g\n", radius);

    if (n == 0) {
        /* Null graph has radius NaN */
        IGRAPH_ASSERT(isnan(radius));
    } else {
        igraph_vector_init(&ecc, 0);
        igraph_eccentricity(graph, weights, &ecc, igraph_vss_vector(center), mode);
        for (igraph_int_t i=0; i < n; i++) {
            IGRAPH_ASSERT(igraph_cmp_epsilon(VECTOR(ecc)[i], radius, eps) == 0);
        }
        igraph_vector_destroy(&ecc);
    }
}

int main(void) {

    igraph_t g;
    igraph_vector_int_t center;
    igraph_vector_t weights, w;

    /* Make sure that this vector has at least as many entries as the
     * largest edge count within the test graphs below. */
    igraph_vector_init_range(&weights, 2, 13);

    igraph_vector_int_init(&center, 0);

    /* Unweighted calculations */
    printf("UNWEIGHTED\n\n");

    printf("Null graph:\n");
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    igraph_graph_center(&g, NULL, &center, IGRAPH_OUT);
    print_vector_int(&center);
    check_radius(&g, &center, IGRAPH_OUT);
    igraph_destroy(&g);

    printf("\nSingleton graph:\n");
    igraph_empty(&g, 1, IGRAPH_UNDIRECTED);
    igraph_graph_center(&g, NULL, &center, IGRAPH_OUT);
    print_vector_int(&center);
    check_radius(&g, &center, IGRAPH_OUT);
    igraph_destroy(&g);

    printf("\nPath with isolated vertex:\n");
    igraph_small(&g, 3, IGRAPH_UNDIRECTED,
                 0,2,
                 -1);
    igraph_graph_center(&g, NULL, &center, IGRAPH_OUT);
    print_vector_int(&center);
    check_radius(&g, &center, IGRAPH_OUT);
    igraph_destroy(&g);

    printf("\nFour isolated vertices:\n");
    igraph_small(&g, 4, IGRAPH_UNDIRECTED, -1);
    igraph_graph_center(&g, NULL, &center, IGRAPH_OUT);
    print_vector_int(&center);
    check_radius(&g, &center, IGRAPH_OUT);
    igraph_destroy(&g);

    printf("\nUndirected path graph P_5:\n");
    igraph_ring(&g, 5, IGRAPH_UNDIRECTED, /* mutual */ false, /* circular */ false);
    igraph_graph_center(&g, NULL, &center, IGRAPH_OUT);
    print_vector_int(&center);
    check_radius(&g, &center, IGRAPH_OUT);
    igraph_destroy(&g);

    printf("\nUndirected graph\n");
    igraph_small(&g, 6, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 0, 2, 3, 0, 4, 1, 2, 5, -1);
    igraph_graph_center(&g, NULL, &center, IGRAPH_OUT);
    print_vector_int(&center);
    check_radius(&g, &center, IGRAPH_OUT);
    igraph_destroy(&g);

    printf("\nDirected path graph P_5:\n");
    igraph_ring(&g, 5, IGRAPH_DIRECTED, /* mutual */ false, /* circular */ false);
    igraph_graph_center(&g, NULL, &center, IGRAPH_OUT);
    print_vector_int(&center);
    check_radius(&g, &center, IGRAPH_OUT);
    igraph_destroy(&g);

    printf("\nUndirected star S_10:\n");
    igraph_star(&g, 10, IGRAPH_STAR_UNDIRECTED, 0);
    igraph_graph_center(&g, NULL, &center, IGRAPH_OUT);
    print_vector_int(&center);
    check_radius(&g, &center, IGRAPH_OUT);
    igraph_destroy(&g);

    printf("\nOut-star S_10:\n");
    igraph_star(&g, 10, IGRAPH_STAR_OUT, 0);
    igraph_graph_center(&g, NULL, &center, IGRAPH_OUT);
    print_vector_int(&center);
    check_radius(&g, &center, IGRAPH_OUT);
    igraph_destroy(&g);

    printf("\nOut-star S_10, undirected mode:\n");
    igraph_star(&g, 10, IGRAPH_STAR_OUT, 0);
    igraph_graph_center(&g, NULL, &center, IGRAPH_ALL);
    print_vector_int(&center);
    check_radius(&g, &center, IGRAPH_ALL);
    igraph_destroy(&g);

    printf("\nIn-star S_10:\n");
    igraph_star(&g, 10, IGRAPH_STAR_IN, 0);
    igraph_graph_center(&g, NULL, &center, IGRAPH_OUT);
    print_vector_int(&center);
    check_radius(&g, &center, IGRAPH_OUT);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    /* Weighted calculations */
    printf("\n\nWEIGHTED\n\n");

    printf("Null graph:\n");
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    w = igraph_vector_view(VECTOR(weights), igraph_ecount(&g));
    igraph_graph_center(&g, &w, &center, IGRAPH_OUT);
    print_vector_int(&center);
    igraph_destroy(&g);

    printf("\nSingleton graph:\n");
    igraph_empty(&g, 1, IGRAPH_UNDIRECTED);
    w = igraph_vector_view(VECTOR(weights), igraph_ecount(&g));
    igraph_graph_center(&g, &w, &center, IGRAPH_OUT);
    print_vector_int(&center);
    check_radius_dijkstra(&g, &w, &center, IGRAPH_OUT);
    igraph_destroy(&g);

    printf("\nPath with isolated vertex:\n");
    igraph_small(&g, 3, IGRAPH_UNDIRECTED,
                 0,2,
                 -1);
    w = igraph_vector_view(VECTOR(weights), igraph_ecount(&g));
    igraph_graph_center(&g, &w, &center, IGRAPH_OUT);
    print_vector_int(&center);
    check_radius_dijkstra(&g, &w, &center, IGRAPH_OUT);
    igraph_destroy(&g);

    printf("\nFour isolated vertices:\n");
    igraph_small(&g, 4, IGRAPH_UNDIRECTED, -1);
    w = igraph_vector_view(VECTOR(weights), igraph_ecount(&g));
    igraph_graph_center(&g, &w, &center, IGRAPH_OUT);
    print_vector_int(&center);
    check_radius_dijkstra(&g, &w, &center, IGRAPH_OUT);
    igraph_destroy(&g);

    printf("\nUndirected path graph P_5:\n");
    igraph_ring(&g, 5, IGRAPH_UNDIRECTED, /* mutual */ false, /* circular */ false);
    w = igraph_vector_view(VECTOR(weights), igraph_ecount(&g));
    igraph_graph_center(&g, &w, &center, IGRAPH_OUT);
    print_vector_int(&center);
    check_radius_dijkstra(&g, &w, &center, IGRAPH_OUT);
    igraph_destroy(&g);

    printf("\nUndirected graph\n");
    igraph_small(&g, 6, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 0, 2, 3, 0, 4, 1, 2, 5, -1);
    w = igraph_vector_view(VECTOR(weights), igraph_ecount(&g));
    igraph_graph_center(&g, &w, &center, IGRAPH_OUT);
    print_vector_int(&center);
    check_radius_dijkstra(&g, &w, &center, IGRAPH_OUT);
    igraph_destroy(&g);

    printf("\nDirected path graph P_5:\n");
    igraph_ring(&g, 5, IGRAPH_DIRECTED, /* mutual */ false, /* circular */ false);
    w = igraph_vector_view(VECTOR(weights), igraph_ecount(&g));
    igraph_graph_center(&g, &w, &center, IGRAPH_OUT);
    print_vector_int(&center);
    check_radius_dijkstra(&g, &w, &center, IGRAPH_OUT);
    igraph_destroy(&g);

    printf("\nUndirected star S_10:\n");
    igraph_star(&g, 10, IGRAPH_STAR_UNDIRECTED, 0);
    w = igraph_vector_view(VECTOR(weights), igraph_ecount(&g));
    igraph_graph_center(&g, &w, &center, IGRAPH_OUT);
    print_vector_int(&center);
    check_radius_dijkstra(&g, &w, &center, IGRAPH_OUT);
    igraph_destroy(&g);

    printf("\nOut-star S_10:\n");
    igraph_star(&g, 10, IGRAPH_STAR_OUT, 0);
    w = igraph_vector_view(VECTOR(weights), igraph_ecount(&g));
    igraph_graph_center(&g, &w, &center, IGRAPH_OUT);
    print_vector_int(&center);
    check_radius_dijkstra(&g, &w, &center, IGRAPH_OUT);
    igraph_destroy(&g);

    printf("\nOut-star S_10, undirected mode:\n");
    igraph_star(&g, 10, IGRAPH_STAR_OUT, 0);
    w = igraph_vector_view(VECTOR(weights), igraph_ecount(&g));
    igraph_graph_center(&g, &w, &center, IGRAPH_ALL);
    print_vector_int(&center);
    check_radius_dijkstra(&g, &w, &center, IGRAPH_ALL);
    igraph_destroy(&g);

    printf("\nIn-star S_10:\n");
    igraph_star(&g, 10, IGRAPH_STAR_OUT, 0);
    w = igraph_vector_view(VECTOR(weights), igraph_ecount(&g));
    igraph_graph_center(&g, &w, &center, IGRAPH_OUT);
    print_vector_int(&center);
    check_radius_dijkstra(&g, &w, &center, IGRAPH_OUT);
    igraph_destroy(&g);

    printf("\nDirected cycle C_5\n");
    igraph_ring(&g, 5, IGRAPH_DIRECTED, /* mutual */ false, /* circular */ true);
    w = igraph_vector_view(VECTOR(weights), igraph_ecount(&g));
    igraph_graph_center(&g, &w, &center, IGRAPH_OUT);
    print_vector_int(&center);
    check_radius_dijkstra(&g, &w, &center, IGRAPH_OUT);

    printf("\nDirected cycle C_5, mode=IN\n");
    igraph_graph_center(&g, &w, &center, IGRAPH_IN);
    print_vector_int(&center);
    check_radius_dijkstra(&g, &w, &center, IGRAPH_IN);
    igraph_destroy(&g);

    igraph_vector_int_destroy(&center);

    igraph_vector_destroy(&weights);

    VERIFY_FINALLY_STACK();

    return 0;
}
