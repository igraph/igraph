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
  igraph_time_t last_v, last_e, last;

  IGRAPH_VECTOR_CONSTANT(edges, 0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,9);
  IGRAPH_VECTOR_TIME_CONSTANT(v_active, 0,0, 1,1, 2,2, 3,3, 4,4);
  IGRAPH_VECTOR_TIME_CONSTANT(e_active, 0, 1, 2, 3, 4, 5, 6, 7, 8);

  IGRAPH_VECTOR_TIME_CONSTANT(v_active2, 0,0, 0,0, 0,0, 0,0, 0,0);
  IGRAPH_VECTOR_TIME_CONSTANT(e_active2, 0, 0, 0, 0, 0, 0, 0, 0, 0);

  igraph_create_temporal(&graph, &edges, 0, IGRAPH_DIRECTED, &e_active, 0,
			 &v_active, 0);

  igraph_time_last(&graph, &last_v, &last_e, &last);
  if (last_v != 4) { return 1; }
  if (last_e != 8) { return 2; }
  if (last   != 8) { return 3; }

  igraph_destroy(&graph);

  /* ------------------------------------ */

  igraph_create_temporal(&graph, &edges, 0, IGRAPH_UNDIRECTED, 0, 0, 0, 0);

  igraph_time_last(&graph, &last_v, &last_e, &last);
  if (last_v != IGRAPH_END) { return 4; }
  if (last_e != IGRAPH_END) { return 5; }
  if (last   != IGRAPH_END) { return 6; }

  igraph_destroy(&graph);

  /* ------------------------------------ */

  igraph_create_temporal(&graph, &edges, 0, IGRAPH_DIRECTED, &e_active2, 0,
			 &v_active2, 0);

  igraph_time_last(&graph, &last_v, &last_e, &last);
  if (last_v != 0) { return 7; }
  if (last_e != 0) { return 8; }
  if (last   != 0) { return 9; }

  igraph_destroy(&graph);

  return 0;
}
