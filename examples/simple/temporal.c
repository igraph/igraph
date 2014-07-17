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

#include <stdlib.h>
#include <igraph.h>

#include "igraph_types_internal.h"


int are_equal(const igraph_t *correct, const igraph_t *graph) {

  if (correct-> n != graph->n) { exit(1); }
  if (correct->directed != graph->directed) { exit(2); }
  if (!igraph_vector_all_e(&correct->from, &graph->from)) { exit(3); }
  if (!igraph_vector_all_e(&correct->to, &graph->to)) { exit(4); }
  if (!igraph_vector_all_e(&correct->oi, &graph->oi)) { exit(5); }
  if (!igraph_vector_all_e(&correct->ii, &graph->ii)) { exit(6); }
  if (!igraph_vector_all_e(&correct->os, &graph->os)) { exit(7); }
  if (!igraph_vector_all_e(&correct->is, &graph->is)) { exit(8); }
  if (!igraph_vector_int_is_null(&graph->vb) &&
      !igraph_vector_int_all_e(&correct->vb, &graph->vb)) { exit(9); }
  if (!igraph_vector_int_is_null(&graph->eb) &&
      !igraph_vector_int_all_e(&correct->eb, &graph->eb)) { exit(10); }
  /* if (!igraph_vector_time_all_e(&correct->vd, &graph->vd)) { exit(11); } */
  /* if (!igraph_vector_time_all_e(&correct->ed, &graph->ed)) { exit(12); } */
  if (correct->now != graph->now) { exit(13); }

  return 0;
}

int check_empty() {
  /* The graph structure */
  igraph_t graph;
  IGRAPH_VECTOR_CONSTANT(edges);
  IGRAPH_VECTOR_TIME_CONSTANT(v_active);
  IGRAPH_VECTOR_TIME_CONSTANT(e_active);

  /* The correct representation */
  IGRAPH_I_STR(correct, /* n= */ 0, IGRAPH_DIRECTED,
	       /* from = */ V(), /* to= */ V(),
	       /* oi = */ V(), /* ii= */ V(),
	       /* os = */ V(0), /* is= */ V(0),
	       /* attr= */ 0,
	       /* vb= */ V(0,0), /* eb= */ V(0,0),
	       /* vd= */ V(), /* ed =*/ V(),
	       /* now= */ -1);

  /* Create graph */
  igraph_create_temporal(&graph, &edges, /* n= */ 0, IGRAPH_DIRECTED,
			 &e_active, /* e_inactive=*/ 0,
			 &v_active, /* v_inactive=*/ 0);

  /* and compare */
  are_equal(&correct, &graph);

  igraph_destroy(&graph);

  return 0;
}

int check_empty5() {

  igraph_t graph;
  IGRAPH_VECTOR_CONSTANT(edges);
  IGRAPH_VECTOR_TIME_CONSTANT(v_active, 0,0,0,0,0);
  IGRAPH_VECTOR_TIME_CONSTANT(e_active);

  IGRAPH_I_STR(correct, 5, IGRAPH_DIRECTED, /* from= */ V(), /* to= */ V(),
	       /* oi= */ V(), /* ii= */ V(), /* os= */ V(0,0,0,0,0,0),
	       /* is= */ V(0,0,0,0,0,0), /* attr= */ 0,
	       /* vb= */ V(0,5), /* eb= */ V(0,0), /* vd= */ V(),
	       /* ed= */ V(), /* now= */ -1);

  igraph_create_temporal(&graph, &edges, igraph_vector_time_size(&v_active),
			 IGRAPH_DIRECTED, &e_active, 0, &v_active, 0);

  are_equal(&correct, &graph);

  igraph_destroy(&graph);

  return 0;
}

int check_empty5_2() {

  igraph_t graph;
  IGRAPH_VECTOR_CONSTANT(edges);
  IGRAPH_VECTOR_TIME_CONSTANT(v_active, 1,1,5,5,7);
  IGRAPH_VECTOR_TIME_CONSTANT(e_active);

  IGRAPH_I_STR(correct, 5, IGRAPH_DIRECTED, /* from= */ V(), /* to= */ V(),
	       /* oi= */ V(), /* ii= */ V(), /* os= */ V(0,0,0,0,0,0),
	       /* is= */ V(0,0,0,0,0,0), /* attr= */ 0,
	       /* vb= */ V(0,0,2,2,2,2,4,4,5), /* eb= */ V(0,0),
	       /* vd= */ V(), /* ed= */ V(), /* now= */ -1);

  igraph_create_temporal(&graph, &edges, igraph_vector_time_size(&v_active),
			 IGRAPH_DIRECTED, &e_active, 0, &v_active, 0);

  are_equal(&correct, &graph);

  igraph_destroy(&graph);

  return 0;
}

int check_edge() {

  igraph_t graph;
  IGRAPH_VECTOR_CONSTANT(edges, 0,1);
  IGRAPH_VECTOR_TIME_CONSTANT(v_active, 0,0,5,5,7);

  IGRAPH_I_STR(correct, 5, IGRAPH_DIRECTED, /* from= */ V(0),
	       /* to= */ V(1), /* oi= */ V(0), /* ii= */ V(0),
	       /* os= */ V(0,1,1,1,1,1), /* is= */ V(0,0,1,1,1,1),
	       /* attr= */ 0, /* vb= */ V(0,2,2,2,2,2,4,4,5),
	       /* eb= */ 0, /* vd= */ V(), /* ed= */ V(),
	       /* now= */ -1);

  igraph_create_temporal(&graph, &edges, igraph_vector_time_size(&v_active),
			 IGRAPH_DIRECTED, 0, 0, &v_active, 0);

  are_equal(&correct, &graph);

  igraph_destroy(&graph);

  return 0;
}

int check_edge_2() {

  igraph_t graph;
  IGRAPH_VECTOR_CONSTANT(edges, 0,1);
  IGRAPH_VECTOR_TIME_CONSTANT(v_active, 0,1,5,5,7);
  IGRAPH_VECTOR_TIME_CONSTANT(e_active, 1);

  IGRAPH_I_STR(correct, 5, IGRAPH_DIRECTED, /* from= */ V(0),
	       /* to= */ V(1), /* oi= */ V(0), /* ii= */ V(0),
	       /* os= */ V(0,1,1,1,1,1), /* is= */ V(0,0,1,1,1,1),
	       /* attr= */ 0, /* vb= */ V(0,1,2,2,2,2,4,4,5),
	       /* eb= */ V(0,0,1), /* vd= */ V(), /* ed= */ V(),
	       /* now= */ -1);

  igraph_create_temporal(&graph, &edges, igraph_vector_time_size(&v_active),
			 IGRAPH_DIRECTED, &e_active, 0, &v_active, 0);

  are_equal(&correct, &graph);

  igraph_destroy(&graph);

  return 0;
}

int check_ring() {

  igraph_t graph;
  IGRAPH_VECTOR_CONSTANT(edges, 0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,0);
  IGRAPH_VECTOR_TIME_CONSTANT(e_active, 9,8,7,6,5,4,3,2,1,0);

  IGRAPH_I_STR(correct, 10, IGRAPH_DIRECTED,
	       /* from= */ V(9,8,7,6,5,4,3,2,1,0),
	       /* to=   */ V(0,9,8,7,6,5,4,3,2,1),
	       /* oi=   */ V(9,8,7,6,5,4,3,2,1,0),
	       /* ii=   */ V(0,9,8,7,6,5,4,3,2,1),
	       /* os=   */ V(0,1,2,3,4,5,6,7,8,9,10),
	       /* is=   */ V(0,1,2,3,4,5,6,7,8,9,10),
	       /* attr= */ 0,
	       /* vb=   */ 0,
	       /* eb=   */ V(0,1,2,3,4,5,6,7,8,9,10),
	       /* vd=   */ V(), /* ed= */ V(), /* now= */ -1);

  igraph_create_temporal(&graph, &edges, 0, IGRAPH_DIRECTED, &e_active, 0,
			 0, 0);

  are_equal(&correct, &graph);

  igraph_destroy(&graph);

  return 0;
}

int check_birth_indexing() {

  igraph_t graph;
  IGRAPH_VECTOR_CONSTANT(edges, 0,1, 0,2, 0,3, 0,4, 0,5,
			 1,6, 1,7, 1,8, 1,9);
  IGRAPH_VECTOR_TIME_CONSTANT(e_active, 2, 5, 9, 1, 4, 3, 6, 8, 7);

  IGRAPH_I_STR(correct, 10, IGRAPH_DIRECTED,
	       /* from= */ V(0,0,1,0,0,1,1,1,0),
	       /* to=   */ V(4,1,6,5,2,7,9,8,3),
	       /* oi=   */ V(0,1,3,4,8,2,5,6,7),
	       /* ii=   */ V(1,4,8,0,3,2,5,7,6),
	       /* os=   */ V(0,5,9,9,9,9,9,9,9,9,9),
	       /* is=   */ V(0,0,1,2,3,4,5,6,7,8,9),
	       /* attr= */ 0,
	       /* vb=   */ 0,
	       /* eb=   */ V(0,0,1,2,3,4,5,6,7,8,9),
	       /* vd=   */ V(),
	       /* ed=   */ V(),
	       /* now=  */ -1);

  igraph_create_temporal(&graph, &edges, 0, IGRAPH_DIRECTED, &e_active, 0,
			 0, 0);

  are_equal(&correct, &graph);

  igraph_destroy(&graph);

  return 0;
}

int main() {
  check_empty();
  check_empty5();
  check_empty5_2();
  check_edge();
  check_edge_2();
  check_ring();
  check_birth_indexing();
  return 0;
}
