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
    IGRAPH_ASSERT(igraph_get_isomorphisms_vf2(g1, g2, vertex_color1, vertex_color2, edge_color1, edge_color2, &maps, node_compat_fn, edge_compat_fn, arg) == error);
    print_and_destroy_maps(&maps);
    printf("\n");
}

void check_print_destroy_simple(igraph_t *g1, igraph_t *g2) {
    check_print_destroy(g1, g2, NULL, NULL, NULL, NULL, NULL, NULL, NULL, IGRAPH_SUCCESS);
}

int main() {
    igraph_t ring, ring_dir;
    igraph_t g_0, g_1;
    igraph_vector_int_t coloring;
    int three = 3;

    igraph_vector_int_init_int(&coloring, 5, 0, 1, 0, 1, 0);
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

    printf("Two rings:\n");
    check_print_destroy_simple(&ring, &ring);

    printf("Two directed rings:\n");
    check_print_destroy_simple(&ring_dir, &ring_dir);

    printf("Two rings where node parity should be equal:\n");
    check_print_destroy(&ring, &ring, NULL, NULL, NULL, NULL, &compat_parity, NULL, NULL, IGRAPH_SUCCESS);

    printf("Two rings where edge parity should be equal:\n");
    check_print_destroy(&ring, &ring, NULL, NULL, NULL, NULL, NULL, &compat_parity, NULL, IGRAPH_SUCCESS);

    printf("Two rings with only one vertex coloring:\n");
    check_print_destroy(&ring, &ring, &coloring, NULL, NULL, NULL, NULL, NULL, NULL, IGRAPH_SUCCESS);

    printf("Two rings with vertex coloring:\n");
    check_print_destroy(&ring, &ring, &coloring, &coloring, NULL, NULL, NULL, NULL, NULL, IGRAPH_SUCCESS);

    printf("Two rings with edge coloring:\n");
    check_print_destroy(&ring, &ring, NULL, NULL, &coloring, &coloring, NULL, NULL, NULL, IGRAPH_SUCCESS);

    printf("Two rings where node of graph 1 should not be 3 higher than node of graph 2:\n");
    check_print_destroy(&ring, &ring, NULL, NULL, NULL, NULL, &compat_not_arg, NULL, &three, IGRAPH_SUCCESS);

    VERIFY_FINALLY_STACK();
    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Two rings with different directedness.\n");
    check_print_destroy(&ring, &ring_dir, NULL, NULL, NULL, NULL, NULL, NULL, NULL, IGRAPH_EINVAL);

    igraph_destroy(&g_0);
    igraph_destroy(&g_1);
    igraph_destroy(&ring);
    igraph_destroy(&ring_dir);
    igraph_vector_int_destroy(&coloring);

    VERIFY_FINALLY_STACK();
    return 0;
}
