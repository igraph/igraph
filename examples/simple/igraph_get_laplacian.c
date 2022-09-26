/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2022  The igraph development team <igraph@igraph.org>

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

int main(void) {
    igraph_t g;
    igraph_vector_t weights;
    igraph_matrix_t m;

    igraph_matrix_init(&m, 1, 1);
    igraph_vector_init_int(&weights, 5, 1, 2, 3, 4, 5);

    igraph_ring(&g, 5, IGRAPH_DIRECTED, 0, 1);
    igraph_get_laplacian(&g, &m, IGRAPH_OUT, IGRAPH_LAPLACIAN_SYMMETRIC, &weights);
    igraph_matrix_print(&m);

    igraph_vector_destroy(&weights);
    igraph_matrix_destroy(&m);
    igraph_destroy(&g);
}
