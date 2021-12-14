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
  igraph_vector_int_t parents;
  igraph_vector_int_init(&parents, 5);

  VECTOR(parents)[0] = 0;
  VECTOR(parents)[1] = 0;
  VECTOR(parents)[2] = 0;
  VECTOR(parents)[3] = 1;
  VECTOR(parents)[4] = 1;

  IGRAPH_ASSERT(igraph_tree_from_parent_vector(&graph, &parents, IGRAPH_TREE_UNDIRECTED) == IGRAPH_SUCCESS);
  igraph_destroy(&graph);

  VERIFY_FINALLY_STACK();
  return 0;
}
