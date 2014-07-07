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

int t1() {

  igraph_t graph;
  igraph_vector_time_t act, inact;

  IGRAPH_VECTOR_CONSTANT(edges, 0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,9);
  IGRAPH_VECTOR_TIME_CONSTANT(v_active, 0,0, 1,1, 2,2, 3,3, 4,4);
  IGRAPH_VECTOR_TIME_CONSTANT(e_active, 0, 1, 2, 3, 4, 5, 6, 7, 8);

  igraph_create_temporal(&graph, &edges, 0, IGRAPH_DIRECTED, &e_active, 0,
			 &v_active, 0);

  igraph_vector_time_init(&act, 0);
  igraph_vector_time_init(&inact, 0);

  igraph_vertices_range(&graph, igraph_vss_all(), &act, &inact);
  if (!igraph_vector_time_all_e(&v_active, &act)) { return 1; }
  if (!igraph_vector_time_isininterval(&inact, IGRAPH_END, IGRAPH_END)) {
    return 2;
  }

  igraph_edges_range(&graph, igraph_ess_all(IGRAPH_EDGEORDER_ID),
		     &act, &inact);
  if (!igraph_vector_time_all_e(&e_active, &act)) { return 3; }
  if (!igraph_vector_time_isininterval(&inact, IGRAPH_END, IGRAPH_END)) {
    return 4;
  }

  igraph_vector_time_destroy(&act);
  igraph_vector_time_destroy(&inact);
  igraph_destroy(&graph);

  return 0;
}

int t2() {

  igraph_t graph;
  igraph_vector_time_t act, inact;
  int i;

  IGRAPH_VECTOR_CONSTANT(edges, 0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,9);
  IGRAPH_VECTOR_TIME_CONSTANT(v_active, 0,0, 1,1, 2,2, 3,3, 4,4);
  IGRAPH_VECTOR_TIME_CONSTANT(e_active, 0, 1, 2, 3, 4, 5, 6, 7, 8);

  igraph_create_temporal(&graph, &edges, 0, IGRAPH_DIRECTED, &e_active, 0,
			 &v_active, 0);

  igraph_vector_time_init(&act, 0);
  igraph_vector_time_init(&inact, 0);

  for (i = 0; i < 10; i++) {
    igraph_vertices_range(&graph, igraph_vss_1(i), &act, &inact);
    if (igraph_vector_time_size(&act) != 1 ||
	VECTOR(act)[0] != VECTOR(v_active)[i]) { return 5; }
    if (igraph_vector_time_size(&inact) != 1 ||
	VECTOR(inact)[0] != IGRAPH_END) { return 6; }
  }

  for (i = 0; i < igraph_ecount(&graph); i++) {
    igraph_edges_range(&graph, igraph_ess_1(i), &act, &inact);
    if (igraph_vector_time_size(&act) != 1 ||
	VECTOR(act)[0] != VECTOR(e_active)[i]) { return 7; }
    if (igraph_vector_time_size(&inact) != 1 ||
	VECTOR(inact)[0] != IGRAPH_END) { return 8; }
  }

  igraph_vector_time_destroy(&act);
  igraph_vector_time_destroy(&inact);
  igraph_destroy(&graph);

  return 0;
}

int main() {
  int ret;
  ret = t1(); if (ret) { return ret; }
  ret = t2(); if (ret) { return ret; }
  return 0;
}
