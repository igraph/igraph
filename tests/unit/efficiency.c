/*
   IGraph library.
   Copyright (C) 2020-2022  The igraph development team <igraph@igraph.org>

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

int test_graph(const char* name, const igraph_t* graph, const igraph_real_t* weights_array) {
    igraph_real_t eff;
    igraph_vector_t eff_vec;
    igraph_vector_t weights;

    printf("###### Testing graph: %s ######\n\n", name);

    igraph_vector_init(&eff_vec, 0);

    if (weights_array) {
        printf("UNWEIGHTED CASE:\n\n");
    }

    igraph_global_efficiency(graph, &eff, NULL, 0);
    printf("Global efficiency, undirected: %f\n", eff);

    igraph_global_efficiency(graph, &eff, NULL, 1);
    printf("Global efficiency, directed: %f\n", eff);

    igraph_average_local_efficiency(graph, &eff, NULL, 0, IGRAPH_ALL);
    printf("Average local efficiency, undirected: %f\n", eff);

    igraph_average_local_efficiency(graph, &eff, NULL, 1, IGRAPH_ALL);
    printf("Average local efficiency, directed, all neighbors: %f\n", eff);

    igraph_average_local_efficiency(graph, &eff, NULL, 1, IGRAPH_IN);
    printf("Average local efficiency, directed, in-neighbors: %f\n", eff);

    igraph_average_local_efficiency(graph, &eff, NULL, 1, IGRAPH_OUT);
    printf("Average local efficiency, directed, out-neighbors: %f\n", eff);

    printf("\nLocal efficiency, undirected:\n");
    igraph_local_efficiency(graph, &eff_vec, igraph_vss_all(), NULL, 0, IGRAPH_ALL);
    print_vector(&eff_vec);

    printf("\nLocal efficiency, directed, all neighbors:\n");
    igraph_local_efficiency(graph, &eff_vec, igraph_vss_all(), NULL, 1, IGRAPH_ALL);
    print_vector(&eff_vec);

    printf("\nLocal efficiency, directed, in-neighbors:\n");
    igraph_local_efficiency(graph, &eff_vec, igraph_vss_all(), NULL, 1, IGRAPH_IN);
    print_vector(&eff_vec);

    printf("\nLocal efficiency, directed, out-neighbors:\n");
    igraph_local_efficiency(graph, &eff_vec, igraph_vss_all(), NULL, 1, IGRAPH_OUT);
    print_vector(&eff_vec);

    if (weights_array) {
        igraph_vector_view(&weights, weights_array, igraph_ecount(graph));
        printf("\nWEIGHTED CASE:\n\n");

        igraph_global_efficiency(graph, &eff, &weights, 0);
        printf("Global efficiency, undirected: %f\n", eff);

        igraph_global_efficiency(graph, &eff, &weights, 1);
        printf("Global efficiency, directed: %f\n", eff);

        igraph_average_local_efficiency(graph, &eff, &weights, 0, IGRAPH_ALL);
        printf("Average local efficiency, undirected: %f\n", eff);

        igraph_average_local_efficiency(graph, &eff, &weights, 1, IGRAPH_ALL);
        printf("Average local efficiency, directed, all neighbors: %f\n", eff);

        igraph_average_local_efficiency(graph, &eff, &weights, 1, IGRAPH_IN);
        printf("Average local efficiency, directed, in-neighbors: %f\n", eff);

        igraph_average_local_efficiency(graph, &eff, &weights, 1, IGRAPH_OUT);
        printf("Average local efficiency, directed, out-neighbors: %f\n", eff);

        printf("\nLocal efficiency, undirected:\n");
        igraph_local_efficiency(graph, &eff_vec, igraph_vss_all(), &weights, 0, IGRAPH_ALL);
        print_vector(&eff_vec);

        printf("\nLocal efficiency, directed, all neighbors:\n");
        igraph_local_efficiency(graph, &eff_vec, igraph_vss_all(), &weights, 1, IGRAPH_ALL);
        print_vector(&eff_vec);

        printf("\nLocal efficiency, directed, in-neighbors:\n");
        igraph_local_efficiency(graph, &eff_vec, igraph_vss_all(), &weights, 1, IGRAPH_IN);
        print_vector(&eff_vec);

        printf("\nLocal efficiency, directed, out-neighbors:\n");
        igraph_local_efficiency(graph, &eff_vec, igraph_vss_all(), &weights, 1, IGRAPH_OUT);
        print_vector(&eff_vec);
    }

    printf("\n\n");

    igraph_vector_destroy(&eff_vec);

    return 0;
}

int test_ring(void) {
    int result;
    igraph_t graph;
    const igraph_real_t weights_array[] = {1, 1, 1, 1};

    igraph_ring(&graph, 4, IGRAPH_DIRECTED, /* mutual = */ 0, /* circular = */ 1);
    result = test_graph("Ring graph", &graph, weights_array);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return result;
}

int test_small_graph(void) {
    int result;
    igraph_t graph;
    const igraph_real_t weights_array[] = {4, 4, 4, 3, 1, 5, 1, 2, 4, 5, 3, 5, 5, 4, 1, 1, 5, 4, 1, 1, 2, 1, 3, 5};

    igraph_small(&graph, 13, /* directed= */ 1,
                 0,8, 1,3, 1,4, 1,5, 1,8, 1,10, 2,0, 2,1, 2,4, 3,5, 4,2, 4,7,
                 4,9, 5,3, 5,10, 6,7, 8,2, 8,3, 8,4, 8,9, 9,3, 9,4, 11,9, 11,3,
                 -1);
    result = test_graph("Small test graph", &graph, weights_array);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return result;
}

int main(void) {

    RUN_TEST(test_ring);
    RUN_TEST(test_small_graph);

    return 0;
}
