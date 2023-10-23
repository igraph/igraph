/*
   IGraph library.
   Copyright (C) 2023  The igraph development team <igraph@igraph.org>

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

int main (void) {
    // TODO: Structure tests and output file
    // Structures that need to be destroyed
    igraph_t g;
    igraph_t g_dir;
    igraph_t g_null;
    igraph_t g_loops;
    igraph_t g_loops_dir;
    igraph_t g_multi;
    igraph_t g_multi_dir;
    igraph_t g_loop_multi;
    igraph_t g_loop_multi_dir;

    igraph_matrix_int_t jdm;
    igraph_matrix_int_t jdm_dir;
    igraph_matrix_int_t jdm_null;
    igraph_matrix_int_t jdm_loops;
    igraph_matrix_int_t jdm_loops_dir;
    igraph_matrix_int_t jdm_multi;
    igraph_matrix_int_t jdm_multi_dir;
    igraph_matrix_int_t jdm_loop_multi;
    igraph_matrix_int_t jdm_loop_multi_dir;

    igraph_vector_int_t weights;

    // Undirected
    igraph_small(&g, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);
    // Directed
    igraph_small(&g_dir, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2, 4, 3,
                 -1);

    // Null graph
    igraph_empty(&g_null, 0, IGRAPH_UNDIRECTED);

    // Undirected with self-loops
    igraph_small(&g_loops, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 1, 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);

    // Directed with self-loops
    igraph_small(&g_loops_dir, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 1, 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2, 4, 3,
                 -1);

    // Multigraph
    igraph_small(&g_multi, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 4, 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);

    // Directed multigraph
    igraph_small(&g_multi_dir, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2, 4, 3, 4, 3,
                 -1);

    // Undirected with self-loops and multiedges
    igraph_small(&g_loop_multi, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 1, 1, 4, 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);

    // Directed with self-loops and multiedges
    igraph_small(&g_loop_multi_dir, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 1, 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2, 4, 3, 4, 3,
                 -1);

    // Initialize a matrix
    igraph_matrix_int_init(&jdm, 1, 1);
    igraph_matrix_int_init(&jdm_dir, 1, 1);
    igraph_matrix_int_init(&jdm_null, 1, 1);
    igraph_matrix_int_init(&jdm_loops, 1, 1);
    igraph_matrix_int_init(&jdm_loops_dir, 1, 1);
    igraph_matrix_int_init(&jdm_multi, 1, 1);
    igraph_matrix_int_init(&jdm_multi_dir, 1, 1);
    igraph_matrix_int_init(&jdm_loop_multi, 1, 1);
    igraph_matrix_int_init(&jdm_loop_multi_dir, 1, 1);

    // Initialize weight vector
    igraph_vector_int_init(&weights, 0);
    igraph_vector_int_print(&weights);

    // Testing types of graphs
    igraph_construct_jdm(&g, &jdm, 4, 4, NULL);
    igraph_construct_jdm(&g_dir, &jdm_dir,-1, 4, NULL);
    igraph_construct_jdm(&g_null, &jdm_null, 4, -1, NULL);
    igraph_construct_jdm(&g_loops, &jdm_loops, -1, 4, NULL);
    igraph_construct_jdm(&g_loops_dir, &jdm_loops_dir, -1, 4, NULL);
    igraph_construct_jdm(&g_multi, &jdm_multi, -1, 4, NULL);
    igraph_construct_jdm(&g_multi_dir, &jdm_multi_dir, -1, 4, NULL);
    igraph_construct_jdm(&g_loop_multi, &jdm_loop_multi, -1, 4, NULL);
    igraph_construct_jdm(&g_loop_multi_dir, &jdm_loop_multi_dir, -1, 4, NULL);

    igraph_matrix_int_print(&jdm);
    igraph_matrix_int_print(&jdm_dir);
    igraph_matrix_int_print(&jdm_null);
    igraph_matrix_int_print(&jdm_loops);
    igraph_matrix_int_print(&jdm_loops_dir);
    igraph_matrix_int_print(&jdm_multi);
    igraph_matrix_int_print(&jdm_multi_dir);
    igraph_matrix_int_print(&jdm_loop_multi);
    igraph_matrix_int_print(&jdm_loop_multi_dir);

    // TODO: Test the resize function
    // -1 in both, -1 in din or dout, 5 in both, 1 in both (should throw error)

    // Clean up
    igraph_destroy(&g);
    igraph_destroy(&g_dir);
    igraph_destroy(&g_null);
    igraph_destroy(&g_loops);
    igraph_destroy(&g_loops_dir);
    igraph_destroy(&g_multi);
    igraph_destroy(&g_multi_dir);
    igraph_destroy(&g_loop_multi);
    igraph_destroy(&g_loop_multi_dir);
    igraph_matrix_int_destroy(&jdm);
    igraph_matrix_int_destroy(&jdm_dir);
    igraph_matrix_int_destroy(&jdm_null);
    igraph_matrix_int_destroy(&jdm_loops);
    igraph_matrix_int_destroy(&jdm_loops_dir);
    igraph_matrix_int_destroy(&jdm_multi);
    igraph_matrix_int_destroy(&jdm_multi_dir);
    igraph_matrix_int_destroy(&jdm_loop_multi);
    igraph_matrix_int_destroy(&jdm_loop_multi_dir);
    igraph_vector_int_destroy(&weights);
    VERIFY_FINALLY_STACK();

    return 0;
}
