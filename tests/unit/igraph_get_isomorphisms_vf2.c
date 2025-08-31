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

/* Vertices/edges with the same parity match */
igraph_bool_t compat_parity(const igraph_t *graph1,
                            const igraph_t *graph2,
                            const igraph_int_t g1_num,
                            const igraph_int_t g2_num,
                            void *arg) {
    IGRAPH_UNUSED(graph1);
    IGRAPH_UNUSED(graph2);
    IGRAPH_UNUSED(arg);
    return (g1_num % 2) == (g2_num % 2);
}

igraph_bool_t compat_not_arg(const igraph_t *graph1,
                             const igraph_t *graph2,
                             const igraph_int_t g1_num,
                             const igraph_int_t g2_num,
                             void *arg) {
    IGRAPH_UNUSED(graph1);
    IGRAPH_UNUSED(graph2);
    IGRAPH_UNUSED(arg);
    return g1_num != *(int*)arg + g2_num;
}

void print_and_destroy_maps(igraph_vector_int_list_t *vp) {
    print_vector_int_list(vp);
    igraph_vector_int_list_destroy(vp);
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
    igraph_vector_int_list_t maps;
    igraph_vector_int_list_init(&maps, 0);
    IGRAPH_ASSERT(igraph_get_isomorphisms_vf2(g1, g2, vertex_color1, vertex_color2, edge_color1, edge_color2, &maps, node_compat_fn, edge_compat_fn, arg) == error);
    print_and_destroy_maps(&maps);
    printf("\n");

}

void check_print_destroy_simple(igraph_t *g1, igraph_t *g2) {
    check_print_destroy(g1, g2, NULL, NULL, NULL, NULL, NULL, NULL, NULL, IGRAPH_SUCCESS);
}

int main(void) {
    igraph_t cycle, cycle_dir, cycle_loop;
    igraph_t g_0, g_1;
    igraph_vector_int_t coloring;
    int three = 3;

    igraph_vector_int_init_int(&coloring, 5, 0, 1, 0, 1, 0);
    igraph_small(&g_0, 0, IGRAPH_UNDIRECTED, -1);
    igraph_small(&g_1, 1, IGRAPH_UNDIRECTED, -1);
    igraph_cycle_graph(&cycle, 5, /*directed*/ false, /*mutual*/ false);
    igraph_cycle_graph(&cycle_dir, 5, /*directed*/ true, /*mutual*/ false);
    igraph_cycle_graph(&cycle_loop, 5, /*directed*/ false, /*mutual*/ false);
    igraph_add_edge(&cycle_loop, 2, 2);

    printf("Two empty graphs:\n");
    check_print_destroy_simple(&g_0, &g_0);

    printf("Two singleton graphs:\n");
    check_print_destroy_simple(&g_1, &g_1);

    printf("Empty and singleton graphs:\n");
    check_print_destroy_simple(&g_0, &g_1);

    printf("Two cycles:\n");
    check_print_destroy_simple(&cycle, &cycle);

    printf("Two directed cycles:\n");
    check_print_destroy_simple(&cycle_dir, &cycle_dir);

    printf("Two cycles where node parity should be equal:\n");
    check_print_destroy(&cycle, &cycle, NULL, NULL, NULL, NULL, &compat_parity, NULL, NULL, IGRAPH_SUCCESS);

    printf("Two cycles where edge parity should be equal:\n");
    check_print_destroy(&cycle, &cycle, NULL, NULL, NULL, NULL, NULL, &compat_parity, NULL, IGRAPH_SUCCESS);

    printf("Two cycles with only one vertex coloring:\n");
    check_print_destroy(&cycle, &cycle, &coloring, NULL, NULL, NULL, NULL, NULL, NULL, IGRAPH_SUCCESS);

    printf("Two cycles with vertex coloring:\n");
    check_print_destroy(&cycle, &cycle, &coloring, &coloring, NULL, NULL, NULL, NULL, NULL, IGRAPH_SUCCESS);

    printf("Two cycles with edge coloring:\n");
    check_print_destroy(&cycle, &cycle, NULL, NULL, &coloring, &coloring, NULL, NULL, NULL, IGRAPH_SUCCESS);

    printf("Two cycles where node of graph 1 should not be 3 higher than node of graph 2:\n");
    check_print_destroy(&cycle, &cycle, NULL, NULL, NULL, NULL, &compat_not_arg, NULL, &three, IGRAPH_SUCCESS);

    VERIFY_FINALLY_STACK();
    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Two cycles with different directedness:\n");
    check_print_destroy(&cycle, &cycle_dir, NULL, NULL, NULL, NULL, NULL, NULL, NULL, IGRAPH_EINVAL);

    printf("Graph with loop edges:\n");
    check_print_destroy(&cycle, &cycle_loop, NULL, NULL, NULL, NULL, NULL, NULL, NULL, IGRAPH_EINVAL);

    igraph_destroy(&g_0);
    igraph_destroy(&g_1);
    igraph_destroy(&cycle);
    igraph_destroy(&cycle_dir);
    igraph_destroy(&cycle_loop);
    igraph_vector_int_destroy(&coloring);

    VERIFY_FINALLY_STACK();
    return 0;
}
