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
  igraph_vector_t neis;
  int i;

  IGRAPH_VECTOR_CONSTANT(edges, 0,1, 0,2, 0,3, 0,4, 0,5,
			 1,6, 1,7, 1,8, 1,9);
  IGRAPH_VECTOR_TIME_CONSTANT(e_active, 2, 5, 9, 1, 4, 3, 6, 8, 7);

  igraph_create_temporal(&graph, &edges, 0, IGRAPH_DIRECTED, &e_active, 0,
			 0, 0);

  igraph_vector_init(&neis, 0);

  /* By default, we are at the end of time */

  igraph_incident(&graph, &neis, 0, IGRAPH_ALL);
  igraph_vector_print(&neis);
  igraph_incident(&graph, &neis, 0, IGRAPH_OUT);
  igraph_vector_print(&neis);

  igraph_incident(&graph, &neis, 1, IGRAPH_ALL);
  igraph_vector_print(&neis);
  igraph_incident(&graph, &neis, 1, IGRAPH_OUT);
  igraph_vector_print(&neis);

  /* Let's rewind the time */

  for (i = 0; i < 10; i++) {
    printf("At %i:\n", i);
    igraph_time_goto(&graph,  i);
    igraph_incident(&graph, &neis, 0, IGRAPH_ALL);
    igraph_vector_print(&neis);
    igraph_incident(&graph, &neis, 0, IGRAPH_OUT);
    igraph_vector_print(&neis);
    igraph_incident(&graph, &neis, 1, IGRAPH_ALL);
    igraph_vector_print(&neis);
    igraph_incident(&graph, &neis, 1, IGRAPH_OUT);
    igraph_vector_print(&neis);
  }

  igraph_vector_destroy(&neis);
  igraph_destroy(&graph);

  /* -------------------------------------------------------------- */

  printf("---\n");

  IGRAPH_VECTOR_CONSTANT(edges2, 1,0, 2,0, 3,0, 4,0, 5,0,
			 6,1, 7,1, 8,1, 9,1);
  IGRAPH_VECTOR_TIME_CONSTANT(e_active2, 2, 5, 9, 1, 4, 3, 6, 8, 7);

  igraph_create_temporal(&graph, &edges2, 0, IGRAPH_DIRECTED, &e_active2,
			 0, 0, 0);

  igraph_vector_init(&neis, 0);

  /* By default, we are at the end of time */

  igraph_incident(&graph, &neis, 0, IGRAPH_ALL);
  igraph_vector_print(&neis);
  igraph_incident(&graph, &neis, 0, IGRAPH_IN);
  igraph_vector_print(&neis);

  igraph_incident(&graph, &neis, 1, IGRAPH_ALL);
  igraph_vector_print(&neis);
  igraph_incident(&graph, &neis, 1, IGRAPH_IN);
  igraph_vector_print(&neis);

  /* Let's rewind the time */

  for (i = 0; i < 10; i++) {
    printf("At %i:\n", i);
    igraph_time_goto(&graph,  i);
    igraph_incident(&graph, &neis, 0, IGRAPH_ALL);
    igraph_vector_print(&neis);
    igraph_incident(&graph, &neis, 0, IGRAPH_IN);
    igraph_vector_print(&neis);
    igraph_incident(&graph, &neis, 1, IGRAPH_ALL);
    igraph_vector_print(&neis);
    igraph_incident(&graph, &neis, 1, IGRAPH_IN);
    igraph_vector_print(&neis);
  }

  igraph_vector_destroy(&neis);
  igraph_destroy(&graph);


  return 0;
}
