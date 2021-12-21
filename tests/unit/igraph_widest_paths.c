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

int main() {

    igraph_integer_t n, m;
    igraph_t g;
    igraph_matrix_t res1, res2;
    igraph_vs_t from, to;
    igraph_vector_t w;

    /* === igraph_widest_paths_dijkstra() === */

    /* 1. Simple Graph */
    n = 5;
    m = 5;
    igraph_small(&g, n, IGRAPH_UNDIRECTED,
                 0, 1, 1, 2, 2, 3, 3, 4, 0, 3,
                 -1);
    igraph_vector_init_real(&w, m, 8.0, 6.0, 10.0, 7.0, 5.0);

    igraph_matrix_init(&res1, 5, 5);
    igraph_matrix_init(&res2, 5, 5);
    igraph_vs_seq(&from, 0, 4);
    igraph_vs_seq(&to, 0, 4);

    igraph_widest_paths_dijkstra(&g, &res1, from, to, &w, IGRAPH_OUT);
    print_matrix_format(&res1, stdout, "%f");

    igraph_widest_paths_floyd_warshall(&g, &res2, from, to, &w, IGRAPH_OUT);
    print_matrix_format(&res2, stdout, "%f");

    IGRAPH_ASSERT(igraph_matrix_all_e(&res1, &res2));

    igraph_vector_destroy(&w);
    igraph_vs_destroy(&to);
    igraph_vs_destroy(&from);
    igraph_matrix_destroy(&res2);
    igraph_matrix_destroy(&res1);
    igraph_destroy(&g);

    /* 2. Graph from Wikipedia */
    n = 7;
    m = 11;
    igraph_small(&g, n, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 1, 2, 1, 3, 2, 4, 2, 5,
                 3, 4, 3, 6, 4, 5, 4, 6, 5, 6,
                 -1);
    igraph_vector_init_real(&w, m, 15.0, 53.0, 40.0, 46.0, 31.0,
                            17.0, 3.0, 11.0, 29.0, 8.0, 40.0);

    igraph_matrix_init(&res1, 1, n);
    igraph_matrix_init(&res2, 1, n);
    igraph_vs_1(&from, 3);
    igraph_vs_seq(&to, 0, n-1);

    igraph_widest_paths_dijkstra(&g, &res1, from, to, &w, IGRAPH_OUT);
    print_matrix_format(&res1, stdout, "%f");

    igraph_widest_paths_floyd_warshall(&g, &res2, from, to, &w, IGRAPH_OUT);
    print_matrix_format(&res2, stdout, "%f");

    IGRAPH_ASSERT(igraph_matrix_all_e(&res1, &res2));

    igraph_vector_destroy(&w);
    igraph_vs_destroy(&to);
    igraph_vs_destroy(&from);
    igraph_matrix_destroy(&res2);
    igraph_matrix_destroy(&res1);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    /* Test Ideas:
    - Duplicates in from / to
    - Directed graphs
    - Behaviour with self loops
    - Behaviour under different modes

    */

    return 0;
}
