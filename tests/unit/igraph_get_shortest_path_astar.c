/* IGraph library.
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
#include "test_utilities.h"

#define LENGTH 10

int check_edges(const igraph_t *graph, const igraph_vector_int_t *vertices,
                const igraph_vector_int_t *edges, int error_code) {

    igraph_bool_t directed = igraph_is_directed(graph);

    igraph_integer_t j, n2 = igraph_vector_int_size(edges);
    if (igraph_vector_int_size(vertices) == 0 && n2 == 0) {
        return 0;
    }
    if (igraph_vector_int_size(vertices) != n2 + 1) {
        exit(error_code + 2);
    }
    for (j = 0; j < n2; j++) {
        igraph_integer_t edge = VECTOR(*edges)[j];
        igraph_integer_t from = VECTOR(*vertices)[j];
        igraph_integer_t to = VECTOR(*vertices)[j + 1];
        if (directed) {
            if (from != IGRAPH_FROM(graph, edge) ||
                to   != IGRAPH_TO  (graph, edge)) {
                exit(error_code);
            }
        } else {
            igraph_integer_t from2 = IGRAPH_FROM(graph, edge);
            igraph_integer_t to2 = IGRAPH_TO(graph, edge);
            igraph_integer_t min1 = from < to ? from : to;
            igraph_integer_t max1 = from < to ? to : from;
            igraph_integer_t min2 = from2 < to2 ? from2 : to2;
            igraph_integer_t max2 = from2 < to2 ? to2 : from2;
            if (min1 != min2 || max1 != max2) {
                exit(error_code + 3);
            }
        }
    }

    return 0;
}

igraph_error_t lattice_heuristic(igraph_real_t *result, igraph_integer_t source_id, igraph_integer_t target_id, void *extra) {
    int x[4];
    for (int i = 0; i < 4; i++) {
        x[i] = source_id % LENGTH;
        source_id /= LENGTH;
    }
    *result = abs(LENGTH - 1 - x[0]) + x[1] + x[2] + x[3];
    return IGRAPH_SUCCESS;
}

struct xyt {
    igraph_vector_t x;
    igraph_vector_t y;
};

igraph_error_t euclidean_heuristic(igraph_real_t *result, igraph_integer_t source_id, igraph_integer_t target_id, void *extra) {
    struct xyt *xyp = extra;
    igraph_real_t xt, xf, yt, yf;
    xt = VECTOR(xyp->x)[target_id];
    yt = VECTOR(xyp->y)[target_id];
    xf = VECTOR(xyp->x)[source_id];
    yf = VECTOR(xyp->y)[source_id];
    *result = sqrt((xt-xf)*(xt-xf) + (yt-yf)*(yt-yf));
    return IGRAPH_SUCCESS;
}

int main(void) {

    igraph_t g;
    igraph_vector_int_t vertices, edges;
    igraph_real_t weights[] = { 1, 2, 3, 4, 5, 1, 1, 1, 1, 1 };
    igraph_vector_t weights_vec;
    struct xyt xy;
    igraph_vector_int_t dimvector;
    igraph_integer_t dims[] = {LENGTH, LENGTH, LENGTH, LENGTH};

    //set seed for grg random graph generation
    igraph_rng_seed(igraph_rng_default(), 42);

    /* Simple ring graph without weights */


    igraph_vector_int_init(&vertices, 0);
    igraph_vector_int_init(&edges, 0);


    printf("Astar on singleton, unweighted, no heuristic:\n");
    igraph_small(&g, 1, IGRAPH_UNDIRECTED, -1);
    igraph_get_shortest_path_astar(&g, /*vertices=*/ &vertices,
                                       /*edges=*/ &edges, /*from=*/ 0, /*to=*/ 0,
                                       /*weights=*/ NULL, /*mode=*/ IGRAPH_OUT,
                                       /*heuristic=*/ NULL, NULL);

    check_edges(&g, &vertices, &edges, /*error code*/10);

    igraph_vector_int_print(&vertices);

    igraph_destroy(&g);

    printf("Astar on ring, unweighted, no heuristic:\n");
    igraph_ring(&g, 10, IGRAPH_UNDIRECTED, 0, 1);
    igraph_get_shortest_path_astar(&g, /*vertices=*/ &vertices,
                                       /*edges=*/ &edges, /*from=*/ 0, /*to=*/ 5,
                                       /*weights=*/ NULL, /*mode=*/ IGRAPH_OUT,
                                       /*heuristic=*/ NULL, NULL);

    check_edges(&g, &vertices, &edges, /*error code*/10);

    igraph_vector_int_print(&vertices);

    /* Same ring, but with weights */

    printf("Astar, weighted, no heuristic:\n");
    igraph_vector_view(&weights_vec, weights, sizeof(weights) / sizeof(weights[0]));
    igraph_get_shortest_path_astar(&g, /*vertices=*/ &vertices,
                                       /*edges=*/ &edges, /*from=*/ 0, /*to=*/ 5,
                                       &weights_vec, IGRAPH_OUT,
                                       /*heuristic=*/ NULL, NULL);

    check_edges(&g, &vertices, &edges, 20);

    igraph_vector_int_print(&vertices);

    igraph_destroy(&g);

    printf("Astar, unweighted, lattice with manhattan distance heuristic:\n");
    igraph_vector_int_view(&dimvector, dims, sizeof(dims)/sizeof(dims[0]));

    igraph_square_lattice(&g, &dimvector, /*nei*/ 1, IGRAPH_UNDIRECTED, /*mutual*/ false, /*periodic*/NULL);

    //print_graph_canon(&g);
    igraph_get_shortest_path_astar(&g, /*vertices=*/ &vertices,
                                       /*edges=*/ &edges, /*from=*/ 0, /*to=*/ LENGTH-1,
                                       /*weights*/NULL, IGRAPH_OUT,
                                       lattice_heuristic, NULL);
    igraph_vector_int_print(&vertices);
    check_edges(&g, &vertices, &edges, /*error code*/60);

    igraph_destroy(&g);

    printf("Astar, unweighted, grg with euclidean distance heuristic:\n");
    igraph_vector_init(&xy.x, 0);
    igraph_vector_init(&xy.y, 0);

    igraph_grg_game(&g, /*nodes*/100, /*radius*/0.2, /*torus*/ false, &xy.x, &xy.y);
    igraph_vector_init(&weights_vec, igraph_ecount(&g));

    for (int i = 0; i < igraph_ecount(&g); i++) {
        igraph_real_t xt, xf, yt, yf;
        xt = VECTOR(xy.x)[IGRAPH_TO(&g, i)];
        xf = VECTOR(xy.x)[IGRAPH_FROM(&g, i)];
        yt = VECTOR(xy.y)[IGRAPH_TO(&g, i)];
        yf = VECTOR(xy.y)[IGRAPH_FROM(&g, i)];
        VECTOR(weights_vec)[i] = sqrt((xt-xf)*(xt-xf) + (yt-yf)*(yt-yf));
    }

    igraph_get_shortest_path_astar(&g, /*vertices=*/ &vertices,
                                       /*edges=*/ &edges, /*from=*/ 0, /*to=*/ LENGTH-1,
                                       /*weights*/&weights_vec, IGRAPH_OUT,
                                       euclidean_heuristic, &xy);

    igraph_vector_int_print(&vertices);
    check_edges(&g, &vertices, &edges, /*error code*/80);

    printf("Check with dijkstra:\n");
    igraph_get_shortest_path_dijkstra(&g, /*vertices=*/ &vertices,
                                       /*edges=*/ &edges, /*from=*/ 0, /*to=*/ LENGTH-1,
                                       /*weights*/&weights_vec, IGRAPH_OUT);

    igraph_vector_int_print(&vertices);

    printf("Same situation but through shortest_path_astar interface:\n");
    igraph_get_shortest_path_astar(&g, &vertices,
                                       &edges, /*from=*/ 0, /*to=*/ LENGTH - 1,
                                       /*weights*/&weights_vec, IGRAPH_OUT,
                                       euclidean_heuristic, &xy);

    igraph_vector_int_print(&vertices);

    printf("Checking from out of range.\n");
    igraph_destroy(&g);
    igraph_small(&g, 2, IGRAPH_UNDIRECTED, 0, 1, -1);
    CHECK_ERROR(igraph_get_shortest_path_astar(&g, NULL, NULL, 10, 1, NULL, IGRAPH_ALL, NULL, NULL), IGRAPH_EINVVID);

    printf("Checking to out of range.\n");
    igraph_destroy(&g);
    igraph_small(&g, 2, IGRAPH_UNDIRECTED, 0, 1, -1);
    CHECK_ERROR(igraph_get_shortest_path_astar(&g, NULL, NULL, 0, 10, NULL, IGRAPH_ALL, NULL, NULL), IGRAPH_EINVVID);

    printf("Checking wrong weight length error.\n");
    igraph_vector_destroy(&weights_vec);
    igraph_vector_init_int(&weights_vec, 0);
    CHECK_ERROR(igraph_get_shortest_path_astar(&g, NULL, NULL, 0, 1, &weights_vec, IGRAPH_ALL, NULL, NULL), IGRAPH_EINVAL);

    printf("Checking negative weight error.\n");
    igraph_vector_destroy(&weights_vec);
    igraph_vector_init_int(&weights_vec, 1, -1);
    CHECK_ERROR(igraph_get_shortest_path_astar(&g, NULL, NULL, 0, 1, &weights_vec, IGRAPH_ALL, NULL, NULL), IGRAPH_EINVAL);

    igraph_vector_int_destroy(&vertices);
    igraph_vector_int_destroy(&edges);
    igraph_vector_int_destroy(&vertices);
    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&weights_vec);
    igraph_vector_destroy(&xy.x);
    igraph_vector_destroy(&xy.y);

    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
