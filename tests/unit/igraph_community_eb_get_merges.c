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

void print_and_destroy(igraph_t *g, igraph_bool_t directed, igraph_vector_int_t *edges,
        igraph_vector_t *weights, igraph_matrix_int_t *res, igraph_vector_int_t *bridges,
        igraph_vector_t *modularity, igraph_vector_int_t *membership) {
    igraph_community_eb_get_merges(g, 1, edges, weights, res, bridges, modularity, membership);
    if (bridges) {
        printf("Bridges:");
        igraph_vector_int_print(bridges);
    }
    if (modularity) {
        printf("Modularity:");
        igraph_vector_print(modularity);
    }
    if (membership) {
        printf("Membership:");
        igraph_vector_int_print(membership);
    }
    if (res) {
        printf("Merges:\n");
        igraph_matrix_int_print(res);
    }
    printf("\n");
    igraph_destroy(g);
}

int main() {
    igraph_t g;
    igraph_vector_int_t edges;
    igraph_vector_t weights;
    igraph_matrix_int_t res;
    igraph_vector_int_t bridges;
    igraph_vector_t modularity;
    igraph_vector_int_t membership;

    igraph_matrix_int_init(&res, 0, 0);
    igraph_vector_init(&weights, 0);
    igraph_vector_int_init(&bridges, 0);
    igraph_vector_init(&modularity, 0);
    igraph_vector_int_init(&membership, 0);

    {
        printf("Graph with no vertices:\n");
        igraph_small(&g, 0, IGRAPH_UNDIRECTED, -1);
        igraph_vector_int_view(&edges, NULL, 0);
        print_and_destroy(&g, 1, &edges, &weights, &res, &bridges, &modularity, &membership);
    }

    {
        printf("Graph with one vertex:\n");
        igraph_small(&g, 0, IGRAPH_UNDIRECTED, -1);
        igraph_vector_int_view(&edges, NULL, 0);
        print_and_destroy(&g, 1, &edges, &weights, &res, &bridges, &modularity, &membership);
    }

    {
        printf("Graph with two vertices, one edge:\n");
        igraph_small(&g, 2, IGRAPH_UNDIRECTED, 0,1, -1);
        igraph_integer_t edge_array[] = {0};
        igraph_vector_int_view(&edges, edge_array, 1);
        print_and_destroy(&g, 1, &edges, NULL, &res, &bridges, &modularity, &membership);
    }

    {
        printf("Triangle, remove one edge:\n");
        igraph_small(&g, 3, IGRAPH_UNDIRECTED, 0,1, 0,2, 1,2, -1);
        igraph_integer_t edge_array[] = {0};
        igraph_vector_int_view(&edges, edge_array, 1);
        print_and_destroy(&g, 1, &edges, NULL, &res, &bridges, &modularity, &membership);
    }
    {
        printf("Triangle, remove one edge, only check merges:\n");
        igraph_small(&g, 3, IGRAPH_UNDIRECTED, 0,1, 0,2, 1,2, -1);
        igraph_integer_t edge_array[] = {0};
        igraph_vector_int_view(&edges, edge_array, 1);
        print_and_destroy(&g, 1, &edges, NULL, &res, NULL, NULL, NULL);
    }

    {
        printf("Triangle, remove two edges:\n");
        igraph_small(&g, 3, IGRAPH_UNDIRECTED, 0,1, 0,2, 1,2, -1);
        igraph_integer_t edge_array[] = {0, 1};
        igraph_vector_int_view(&edges, edge_array, 2);
        print_and_destroy(&g, 1, &edges, NULL, &res, &bridges, &modularity, &membership);
    }

    igraph_matrix_int_destroy(&res);
    igraph_vector_destroy(&weights);
    igraph_vector_int_destroy(&bridges);
    igraph_vector_destroy(&modularity);
    igraph_vector_int_destroy(&membership);

    VERIFY_FINALLY_STACK();
}
