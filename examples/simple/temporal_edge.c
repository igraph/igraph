/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2014  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

int main() {

  igraph_t graph;
  igraph_integer_t from, to;
  int ret;

  IGRAPH_VECTOR_CONSTANT(edges, 0,1, 0,2, 0,3, 0,4, 0,5,
			 1,6, 1,7, 1,8, 1,9);
  IGRAPH_VECTOR_TIME_CONSTANT(e_active, 2, 5, 9, 1, 4, 3, 6, 8, 7);

  igraph_create_temporal(&graph, &edges, 0, IGRAPH_DIRECTED, &e_active, 0,
			 0, 0);

  igraph_set_error_handler(igraph_error_handler_ignore);

  igraph_time_goto(&graph, 0);
  ret = igraph_edge(&graph, 0, &from, &to);
  if (ret != IGRAPH_EINVAL) { return 1; }

  igraph_time_goto(&graph, 1);
  ret = igraph_edge(&graph, 0, &from, &to);
  if (ret != IGRAPH_SUCCESS) { return 2; }

  ret = igraph_edge(&graph, 1, &from, &to);
  if (ret != IGRAPH_EINVAL) { return 3; }

  igraph_destroy(&graph);

  return 0;
}
