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
  igraph_vector_time_t active;

  IGRAPH_VECTOR_CONSTANT(edges, 0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,9);
  IGRAPH_VECTOR_TIME_CONSTANT(v_active, 0,0, 1,1, 2,2, 3,3, 4,4);
  IGRAPH_VECTOR_TIME_CONSTANT(e_active, 0, 1, 2, 3, 4, 5, 6, 7, 8);
  IGRAPH_VECTOR_CONSTANT(new_edges, 1,0, 1,0);
  IGRAPH_VECTOR_TIME_CONSTANT(e_active_new, 3, 1);

  igraph_create_temporal(&graph, &edges, 0, IGRAPH_DIRECTED, &e_active, 0,
			 &v_active, 0);

  igraph_add_edges_at(&graph, &new_edges, &e_active_new, 0, 0);

  igraph_add_edges_at(&graph, &new_edges, 0, 0, 0);

  igraph_time_goto(&graph, 2);
  igraph_add_edges_at(&graph, &new_edges, 0, 0, 0);

  igraph_time_goto(&graph, 0);
  igraph_add_edges(&graph, &new_edges, 0);

  igraph_vector_time_init(&active, 0);
  igraph_edges_range(&graph, igraph_ess_all(IGRAPH_EDGEORDER_ID),
		     &active, 0);
  igraph_vector_time_print(&active);

  igraph_time_goto(&graph, IGRAPH_END);

  igraph_vector_time_destroy(&active);
  igraph_destroy(&graph);

  return 0;
}
