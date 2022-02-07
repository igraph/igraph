/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2022 The igraph development team

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

#include "test_utilities.inc"

void test_graphlets_candidate_basis_simple() {
    igraph_t g;
    igraph_vector_int_list_t cliques;
    igraph_vector_t thresholds;

    igraph_full(&g, 5, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_vector_int_list_init(&cliques, 0);
    igraph_vector_init(&thresholds, 0);

    igraph_graphlets_candidate_basis(&g, 0, &cliques, &thresholds);

    igraph_vector_int_list_destroy(&thresholds);
    igraph_vector_destroy(&thresholds);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
}

int main() {
    test_graphlet_basis_simple();

    return 0;
}
