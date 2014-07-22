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
  igraph_vector_t vec;

  IGRAPH_VECTOR_CONSTANT(edges, 0,1, 0,2, 0,3, 0,4, 0,5,
			 1,6, 1,7, 1,8, 1,9);
  IGRAPH_VECTOR_TIME_CONSTANT(v_active, 0, 2, 5, 9, 0, 1, 0, 0, 8, 7);
  IGRAPH_VECTOR_TIME_CONSTANT(e_active, 2, 5, 9, 1, 4, 3, 6, 8, 7);
  IGRAPH_VECTOR_CONSTANT(birth, 0, 0, 0, 0, 1, 2, 5, 7, 8, 9);
  IGRAPH_VECTOR_TIME_CONSTANT(v2_active, 0, 5, 10);

  igraph_i_set_attribute_table(&igraph_cattribute_table);

  igraph_create_temporal(&graph, &edges, 0, IGRAPH_DIRECTED,
			 &e_active, 0, &v_active, 0);

  SETGAN(&graph, "id", 10);
  SETGAS(&graph, "name", "foobar");
  SETVANV(&graph, "name", &birth);

  igraph_time_goto(&graph, 3);

  igraph_add_vertices_at(&graph, 3, &v2_active, /* v_inactive= */ 0,
			 /* attr= */ 0);

  igraph_vector_init(&vec, 0);
  VANV(&graph, "name", &vec);
  igraph_vector_print(&vec);

  igraph_destroy(&graph);

  return 0;
}
