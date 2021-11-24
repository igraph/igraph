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

    igraph_t graph;
    igraph_vector_int_t l;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&l, 0);
    igraph_vector_int_push_back(&l, 1);
    igraph_vector_int_push_back(&l, 2);

    IGRAPH_ASSERT(igraph_circulant(&graph, /* n */ 5, /* l */ &l, /* directed */ 1) == IGRAPH_SUCCESS);

    // todo: add more tests

    return 0;
}