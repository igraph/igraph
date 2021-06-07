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
#include "test_utilities.inc"

void check(igraph_t *graph, igraph_vs_t *vs) {
    igraph_vit_t vit;

    IGRAPH_ASSERT(igraph_vit_create(graph, *vs, &vit) == IGRAPH_SUCCESS);
    for (; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
        printf("%" IGRAPH_PRId "\n", IGRAPH_VIT_GET(vit));
    }
}

int main() {
    igraph_t g, g_no_vertices, g_no_edges;
    igraph_vs_t vs;
    igraph_vector_t v;
    igraph_vit_t vit;

    igraph_small(&g, 5, IGRAPH_DIRECTED, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, -1);
    igraph_small(&g_no_vertices, 0, IGRAPH_UNDIRECTED, -1);
    igraph_small(&g_no_edges, 5, IGRAPH_UNDIRECTED, -1);

    printf("Checking vs_none vertex selector:\n");
    IGRAPH_ASSERT(igraph_vs_none(&vs) == IGRAPH_SUCCESS);
    check(&g, &vs);
    check(&g_no_edges, &vs);
    check(&g_no_vertices, &vs);

    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Checking vector selector:\n");
    igraph_vector_init_int(&v, 3, 2, 3, 4);
    IGRAPH_ASSERT(igraph_vs_vector(&vs, &v) == IGRAPH_SUCCESS);
    printf("Some graph:\n");
    check(&g, &vs);
    printf("Edgeless graph:\n");
    check(&g_no_edges, &vs);
    printf("Graph without vertices should fail.\n");
    IGRAPH_ASSERT(igraph_vit_create(&g_no_vertices, vs, &vit) == IGRAPH_EINVVID);
    igraph_vector_destroy(&v);

    printf("Vertex selector with negative index should fail\n");
    igraph_vector_init_int(&v, 3, -2, 3, 4);
    IGRAPH_ASSERT(igraph_vs_vector(&vs, &v) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vit_create(&g, vs, &vit) == IGRAPH_EINVVID);
    igraph_vector_destroy(&v);

    printf("Checking copy vector selector:\n");
    igraph_vector_init_int(&v, 3, 2, 3, 4);
    IGRAPH_ASSERT(igraph_vs_vector_copy(&vs, &v) == IGRAPH_SUCCESS);
    printf("Some graph:\n");
    check(&g, &vs);
    printf("Edgeless graph:\n");
    check(&g_no_edges, &vs);
    printf("Graph without vertices should fail.\n");
    IGRAPH_ASSERT(igraph_vit_create(&g_no_vertices, vs, &vit) == IGRAPH_EINVVID);
    IGRAPH_ASSERT(igraph_vs_type(&vs) == IGRAPH_VS_VECTOR);
    igraph_vector_destroy(&v);

    igraph_destroy(&g);
    igraph_destroy(&g_no_vertices);
    igraph_destroy(&g_no_edges);

    VERIFY_FINALLY_STACK();
    return 0;
}
