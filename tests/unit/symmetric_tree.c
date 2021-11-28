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

    igraph_t g;
    igraph_vector_int_t v;
    igraph_vector_int_init(&v, 3);
    VECTOR(v)[0] = 3;
    VECTOR(v)[1] = 4;
    VECTOR(v)[2] = 5;


    IGRAPH_ASSERT(igraph_symmetric_tree(&g, &v, IGRAPH_TREE_OUT) == IGRAPH_SUCCESS);

    igraph_vector_int_init(&v, 4);
    VECTOR(v)[0] = 3;
    VECTOR(v)[1] = 4;
    VECTOR(v)[2] = 5;
    VECTOR(v)[3] = 6;

    IGRAPH_ASSERT(igraph_symmetric_tree(&g, &v, IGRAPH_TREE_OUT) == IGRAPH_FAILURE);
    // TODO: add more tests

     VERIFY_FINALLY_STACK();
    return 0;
}
