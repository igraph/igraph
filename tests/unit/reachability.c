/*
   igraph library.
   Copyright (C) 2024  The igraph development team <igraph@igraph.org>

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

void compute_and_print(const igraph_t *g, igraph_neimode_t mode) {
    igraph_vector_int_t membership, csize, reach_counts;
    igraph_bitset_list_t reach;
    igraph_int_t no_of_nodes, no_of_components;

    no_of_nodes = igraph_vcount(g);

    igraph_vector_int_init(&membership, 0);
    igraph_vector_int_init(&csize, 0);
    igraph_bitset_list_init(&reach, 0);

    igraph_reachability(g, &membership, &csize, &no_of_components, &reach, mode);

    igraph_vector_int_init(&reach_counts, no_of_nodes);

    igraph_count_reachable(g, &reach_counts, mode);

    printf("Mode: ");
    switch (mode) {
    case IGRAPH_OUT:
        printf("OUT\n"); break;
    case IGRAPH_IN:
        printf("IN\n"); break;
    case IGRAPH_ALL:
        printf("ALL\n"); break;
    }
    print_vector_int(&membership);
    print_vector_int(&csize);
    printf("No. of components: %" IGRAPH_PRId "\n", no_of_components);
    print_bitset_list(&reach);
    print_vector_int(&reach_counts);

    igraph_bitset_list_destroy(&reach);
    igraph_vector_int_destroy(&csize);
    igraph_vector_int_destroy(&membership);
    igraph_vector_int_destroy(&reach_counts);
}

int main(void) {
    igraph_t g;

    /* Component calculations are run twice to exercise the cache */

    printf("\nNull graph (not connected)\n");
    igraph_empty(&g, 0, IGRAPH_DIRECTED);
    compute_and_print(&g, IGRAPH_OUT);
    compute_and_print(&g, IGRAPH_ALL);
    igraph_destroy(&g);

    printf("\nSingleton graph (connected)\n");
    igraph_empty(&g, 1, IGRAPH_DIRECTED);
    compute_and_print(&g, IGRAPH_OUT);
    compute_and_print(&g, IGRAPH_ALL);
    igraph_destroy(&g);

    printf("\nKautz graph (connected)\n");
    igraph_kautz(&g, 2, 2);
    compute_and_print(&g, IGRAPH_OUT);
    compute_and_print(&g, IGRAPH_ALL);
    igraph_destroy(&g);

    printf("\nDirected 2-path\n");
    igraph_small(&g, 2, IGRAPH_DIRECTED, 0, 1, -1);
    compute_and_print(&g, IGRAPH_OUT);
    compute_and_print(&g, IGRAPH_ALL);
    igraph_destroy(&g);

    printf("\nTwo disjoint 3-cycles\n");
    igraph_small(&g, 6, IGRAPH_DIRECTED, 0, 1, 1, 2, 2, 0, 3, 4, 4, 5, 5, 3, -1);
    compute_and_print(&g, IGRAPH_OUT);
    compute_and_print(&g, IGRAPH_IN);
    compute_and_print(&g, IGRAPH_ALL);
    igraph_destroy(&g);

    printf("\nPath graph with 6 vertices ascending\n");
    igraph_small(&g, 6, IGRAPH_DIRECTED, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, -1);
    compute_and_print(&g, IGRAPH_OUT);
    compute_and_print(&g, IGRAPH_IN);
    compute_and_print(&g, IGRAPH_ALL);
    igraph_destroy(&g);

    printf("\nPath graph with 6 vertices descending\n");
    igraph_small(&g, 6, IGRAPH_DIRECTED, 5, 4, 4, 3, 3, 2, 2, 1, 1, 0, -1);
    compute_and_print(&g, IGRAPH_OUT);
    compute_and_print(&g, IGRAPH_IN);
    compute_and_print(&g, IGRAPH_ALL);
    igraph_destroy(&g);

    printf("\nSmall directed graph\n");
    igraph_small(&g, 13, IGRAPH_DIRECTED, 0, 1, 1, 2, 2, 0, 1, 3, 3, 4, 4, 5, 5, 4, 7, 4, 7, 8, 9, 8, 10, 9, 8, 10, 11, 6, 12, 6, -1);
    compute_and_print(&g, IGRAPH_OUT);
    compute_and_print(&g, IGRAPH_IN);
    compute_and_print(&g, IGRAPH_ALL);
    igraph_destroy(&g);

    printf("\nSmall undirected graph\n");
    igraph_small(&g, 6, IGRAPH_UNDIRECTED,
                 0,1, 1,2, 0,2,
                 3,4,
                 -1);
    compute_and_print(&g, IGRAPH_ALL);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
