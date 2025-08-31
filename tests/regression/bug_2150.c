/*
   igraph library.
   Copyright (C) 2022  The igraph development team

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

#include "../unit/test_utilities.h"

/* Regression test for https://github.com/igraph/igraph/issues/2150 */

int main(void) {
    igraph_t graph;
    igraph_vector_t w;
    igraph_vector_int_list_t res;

    igraph_full(&graph, 3, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_vector_init_int(&w, 3, 1, 2, 3);
    igraph_vector_int_list_init(&res, 0);

    igraph_weighted_cliques(&graph, &w, &res, false, 1, 3, IGRAPH_UNLIMITED);

    igraph_vector_int_list_destroy(&res);
    igraph_vector_destroy(&w);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
