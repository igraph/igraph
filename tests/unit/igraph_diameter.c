/* -*- mode: C -*-  */
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

int main(void) {

    igraph_t g;
    igraph_real_t result;
    igraph_integer_t from, to;
    igraph_vector_int_t path_edge, path_vertex, edge_vec;
    igraph_vector_t weights_vec;
    igraph_real_t weights[] = { 1, 2, 3, 4, 5, 1, 1, 1, 1};
    igraph_integer_t vec[] = {2,8};
    igraph_es_t edge_sele;

    printf("diameter of Barabasi graph:\n");

    igraph_rng_seed(igraph_rng_default(), 1234);
    igraph_barabasi_game(&g, 30, /*power=*/ 1, 30, 0, 0, /*A=*/ 1,
                         IGRAPH_DIRECTED, IGRAPH_BARABASI_BAG,
                         /*start_from=*/ 0);
    igraph_diameter(&g, &result, NULL, NULL, NULL, NULL, IGRAPH_UNDIRECTED, 1);

    printf("Diameter: %" IGRAPH_PRId "\n", (igraph_integer_t) result);

    igraph_destroy(&g);

    printf("diameter of ring and the path in terms of edges and vertices \n");

    igraph_ring(&g, 10, IGRAPH_DIRECTED, 0, 0);
    igraph_vector_int_init(&path_vertex, 0);
    igraph_vector_int_init(&path_edge, 0);
    igraph_vector_view(&weights_vec, weights, sizeof(weights) / sizeof(weights[0]));
    igraph_diameter(&g, &result, &from, &to, &path_vertex, &path_edge, IGRAPH_DIRECTED, 1);
    printf("diameter: %g, from %" IGRAPH_PRId " to %" IGRAPH_PRId "\n", result,
            from, to);
    print_vector_int(&path_vertex);
    print_vector_int(&path_edge);
    igraph_vector_int_destroy(&path_vertex);
    igraph_vector_int_destroy(&path_edge);

    //disconnected graph
    printf("disconnected ring graph\n");
    igraph_vector_int_view(&edge_vec, vec, sizeof(vec) / sizeof(vec[0]));
    igraph_es_vector(&edge_sele, &edge_vec);
    igraph_delete_edges(&g, edge_sele);
    printf("The largest path in one connected component\n");
    igraph_diameter(&g, &result, NULL, NULL, NULL, NULL, IGRAPH_DIRECTED, 1);
    print_real(stdout, result, "%g");
    printf("\nuconn = False \n");
    igraph_diameter(&g, &result, NULL, NULL, NULL, NULL, IGRAPH_DIRECTED, 0);
    print_real(stdout, result, "%g");

    igraph_es_destroy(&edge_sele);
    igraph_destroy(&g);

    // test graph with zero nodes
    printf("\ngraph with zero nodes\n");
    igraph_empty(&g, 0, IGRAPH_DIRECTED);
    igraph_diameter(&g, &result, &from, &to, NULL, NULL, IGRAPH_DIRECTED, 1);
    print_real(stdout, result, "%g");
    printf("\nfrom = %" IGRAPH_PRId ", to = %" IGRAPH_PRId "\n", from, to);
    igraph_destroy(&g);

    //test graph with one node
    printf("graph with one node\n");
    igraph_empty(&g, 1, IGRAPH_DIRECTED);
    igraph_diameter(&g, &result, &from, &to, NULL, NULL, IGRAPH_DIRECTED, 1);
    print_real(stdout, result, "%g");
    printf("\nfrom = %" IGRAPH_PRId ", to = %" IGRAPH_PRId "\n", from, to);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
