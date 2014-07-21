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
  igraph_vector_t edges2;
  igraph_integer_t i, eid, ec;

  IGRAPH_VECTOR_CONSTANT(edges, 0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,9);
  IGRAPH_VECTOR_TIME_CONSTANT(e_active, 2, 4, 5, 8, 6, 3, 1, 0, 7);

  igraph_create_temporal(&graph, &edges, 0, IGRAPH_DIRECTED, &e_active, 0,
			 0, 0);

  igraph_vector_init(&edges2, 0);
  igraph_get_edgelist(&graph, &edges2, /* bycol= */ 0);
  igraph_vector_print(&edges2);

  ec = igraph_ecount(&graph);

  for (i = 0; i < ec; i++) {
    int from = VECTOR(edges2)[ 2 * i ];
    int to   = VECTOR(edges2)[ 2 * i + 1];
    igraph_get_eid(&graph, &eid, from, to,
		   /* directed= */ 0, /* error= */ 1);
    if (eid != i) { return i + 1; }
  }

  igraph_time_goto(&graph, 4);
  for (i = 0 ; i <= 4; i++) {
    int from = VECTOR(edges2)[ 2 * i ];
    int to   = VECTOR(edges2)[ 2 * i + 1];
    igraph_get_eid(&graph, &eid, from, to,
		   /* directed= */ 0, /* error= */ 1);
    if (eid != i) { return ec + i + 1; }
  }

  for (i = 5 ; i < ec; i++) {
    int from = VECTOR(edges2)[ 2 * i ];
    int to   = VECTOR(edges2)[ 2 * i + 1];
    igraph_get_eid(&graph, &eid, from, to,
		   /* directed= */ 0, /* error= */ 0);
    if (eid != -1) { return ec + i + 1; }
  }

  igraph_time_goto(&graph, IGRAPH_END);
  igraph_get_eids(&graph, &edges2, &edges, /* path= */ 0,
		  /* directed= */ 1, /* error= */ 1);
  igraph_vector_print(&edges2);

  igraph_time_goto(&graph, 4);
  igraph_get_eids(&graph, &edges2, &edges, /* path= */ 0,
		  /* directed= */ 1, /* error= */ 0);
  igraph_vector_print(&edges2);

  igraph_time_goto(&graph, IGRAPH_END);
  igraph_get_eids(&graph, &edges2, &edges, /* path= */ 0,
		  /* directed= */ 0, /* error= */ 1);
  igraph_vector_print(&edges2);

  igraph_time_goto(&graph, 4);
  igraph_get_eids(&graph, &edges2, &edges, /* path= */ 0,
		  /* directed= */ 0, /* error= */ 0);
  igraph_vector_print(&edges2);

  igraph_destroy(&graph);

  return 0;
}
