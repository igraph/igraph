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

/* Vertices/edges with the same parity match */
igraph_bool_t compat_parity(const igraph_t *graph1,
                            const igraph_t *graph2,
                            const igraph_integer_t g1_num,
                            const igraph_integer_t g2_num,
                            void *arg) {
    IGRAPH_UNUSED(graph1);
    IGRAPH_UNUSED(graph2);
    IGRAPH_UNUSED(arg);
    return (g1_num % 2) == (g2_num % 2);
}

igraph_bool_t compat_not_arg(const igraph_t *graph1,
                             const igraph_t *graph2,
                             const igraph_integer_t g1_num,
                             const igraph_integer_t g2_num,
                             void *arg) {
    IGRAPH_UNUSED(graph1);
    IGRAPH_UNUSED(graph2);
    IGRAPH_UNUSED(arg);
    return g1_num != *(int*)arg + g2_num;
}

void print_and_destroy_maps(igraph_vector_ptr_t *vp) {
    long int i;
    for (i = 0; i < igraph_vector_ptr_size(vp); i++) {
        print_vector(VECTOR(*vp)[i]);
        igraph_vector_destroy(VECTOR(*vp)[i]);
        igraph_free(VECTOR(*vp)[i]);
    }
    igraph_vector_ptr_destroy(vp);
}

void check_print_destroy(igraph_t *g1,
                         igraph_t *g2,
                         igraph_vector_int_t *vertex_color1,
                         igraph_vector_int_t *vertex_color2,
                         igraph_vector_int_t *edge_color1,
                         igraph_vector_int_t *edge_color2,
                         igraph_isocompat_t *node_compat_fn,
                         igraph_isocompat_t *edge_compat_fn,
                         void *arg,
                         int error) {
    igraph_vector_ptr_t maps;
    igraph_vector_ptr_init(&maps, 0);
    IGRAPH_ASSERT(igraph_get_subisomorphisms_vf2(g1, g2, vertex_color1, vertex_color2, edge_color1, edge_color2, &maps, node_compat_fn, edge_compat_fn, arg) == error);
    print_and_destroy_maps(&maps);
    printf("\n");
}

void check_print_destroy_simple(igraph_t *g1, igraph_t *g2) {
    check_print_destroy(g1, g2, NULL, NULL, NULL, NULL, NULL, NULL, NULL, IGRAPH_SUCCESS);
}

int main() {
    igraph_t ring, ring_dir;
    igraph_t ring_plus, ring_plus_dir;
    igraph_t g_0, g_1;
    igraph_vector_int_t coloring;
    igraph_vector_int_t plus_edge_coloring;
    igraph_vector_int_t plus_vertex_coloring;
    int three = 3;

    igraph_vector_int_init_int(&plus_vertex_coloring, 6, 0, 1, 0, 1, 1, 0);
    igraph_vector_int_init_int(&coloring, 5, 0, 1, 0, 1, 0);
    igraph_vector_int_init_int(&plus_edge_coloring, 6, 0, 0, 1, 0, 1, 0);
    igraph_small(&ring_plus, 6, 0, 0,1, 0,3, 2,1, 2,3, 3,5, 5,0, -1);
    igraph_small(&ring_plus_dir, 6, 1, 0,1, 0,3, 1,2, 2,3, 3,5, 5,0, -1);
    igraph_small(&g_0, 0, 0, -1);
    igraph_small(&g_1, 1, 0, -1);
    igraph_ring(&ring, 5, /*directed*/ 0, /*mutual*/ 0, /*circular*/ 1);
    igraph_ring(&ring_dir, 5, /*directed*/ 1, /*mutual*/ 0, /*circular*/ 1);

    printf("Two empty graphs:\n");
    check_print_destroy_simple(&g_0, &g_0);

    printf("Two singleton graphs:\n");
    check_print_destroy_simple(&g_1, &g_1);

    printf("Empty and singleton graphs:\n");
    check_print_destroy_simple(&g_0, &g_1);

    printf("Singleton and empty graphs:\n");
    check_print_destroy_simple(&g_1, &g_0);

    printf("Ring with add vertex and edge (ring+) and ring:\n");
    check_print_destroy_simple(&ring_plus, &ring);

    printf("Ring+ and ring, directed:\n");
    check_print_destroy_simple(&ring_plus_dir, &ring_dir);

    printf("Ring+ and ring where node parity should be equal:\n");
    check_print_destroy(&ring_plus, &ring, NULL, NULL, NULL, NULL, &compat_parity, NULL, NULL, IGRAPH_SUCCESS);

    printf("Ring+ and ring where edge parity should be equal:\n");
    check_print_destroy(&ring_plus, &ring, NULL, NULL, NULL, NULL, NULL, &compat_parity, NULL, IGRAPH_SUCCESS);

    printf("Ring+ and ring with only one vertex coloring:\n");
    check_print_destroy(&ring_plus, &ring, &coloring, NULL, NULL, NULL, NULL, NULL, NULL, IGRAPH_SUCCESS);

    printf("Ring+ and ring with vertex coloring:\n");
    check_print_destroy(&ring_plus, &ring, &plus_vertex_coloring, &coloring, NULL, NULL, NULL, NULL, NULL, IGRAPH_SUCCESS);

    printf("Ring+ and ring with edge coloring:\n");
    check_print_destroy(&ring_plus, &ring, NULL, NULL, &plus_edge_coloring, &coloring, NULL, NULL, NULL, IGRAPH_SUCCESS);

    printf("Ring+ and ring where node of graph 1 should not be 3 higher than node of graph 2:\n");
    check_print_destroy(&ring_plus, &ring, NULL, NULL, NULL, NULL, &compat_not_arg, NULL, &three, IGRAPH_SUCCESS);

    VERIFY_FINALLY_STACK();
    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Ring+ and ring with different directedness.\n");
    check_print_destroy(&ring_plus_dir, &ring, NULL, NULL, NULL, NULL, NULL, NULL, NULL, IGRAPH_EINVAL);

    igraph_destroy(&g_0);
    igraph_destroy(&g_1);
    igraph_destroy(&ring);
    igraph_destroy(&ring_dir);
    igraph_destroy(&ring_plus);
    igraph_destroy(&ring_plus_dir);
    igraph_vector_int_destroy(&coloring);
    igraph_vector_int_destroy(&plus_edge_coloring);
    igraph_vector_int_destroy(&plus_vertex_coloring);

    VERIFY_FINALLY_STACK();
    return 0;
}
