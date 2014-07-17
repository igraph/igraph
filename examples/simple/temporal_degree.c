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
  igraph_vector_t deg;
  int i;

  IGRAPH_VECTOR_CONSTANT(edges, 0,1, 0,2, 0,3, 0,4, 0,5,
			 1,6, 1,7, 1,8, 1,9);
  IGRAPH_VECTOR_TIME_CONSTANT(e_active, 2, 5, 9, 1, 4, 3, 6, 8, 7);

  igraph_create_temporal(&graph, &edges, 0, IGRAPH_DIRECTED, &e_active, 0,
			 0, 0);

  igraph_vector_init(&deg, 0);

  printf("ALL, LOOPS:\n");
  for (i = 0; i < 10; i++) {
    igraph_time_goto(&graph, i);
    igraph_degree(&graph, &deg, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    printf("  %i: ", i); igraph_vector_print(&deg);
  }

  printf("ALL, NO_LOOPS:\n");
  for (i = 0; i < 10; i++) {
    igraph_time_goto(&graph, i);
    igraph_degree(&graph, &deg, igraph_vss_all(), IGRAPH_ALL,
		  IGRAPH_NO_LOOPS);
    printf("  %i: ", i); igraph_vector_print(&deg);
  }

  printf("OUT, LOOPS:\n");
  for (i = 0; i < 10; i++) {
    igraph_time_goto(&graph, i);
    igraph_degree(&graph, &deg, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS);
    printf("  %i: ", i); igraph_vector_print(&deg);
  }

  printf("OUT, NO_LOOPS:\n");
  for (i = 0; i < 10; i++) {
    igraph_time_goto(&graph, i);
    igraph_degree(&graph, &deg, igraph_vss_all(), IGRAPH_OUT,
		  IGRAPH_NO_LOOPS);
    printf("  %i: ", i); igraph_vector_print(&deg);
  }

  printf("IN, LOOPS:\n");
  for (i = 0; i < 10; i++) {
    igraph_time_goto(&graph, i);
    igraph_degree(&graph, &deg, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS);
    printf("  %i: ", i); igraph_vector_print(&deg);
  }

  printf("IN, NO_LOOPS:\n");
  for (i = 0; i < 10; i++) {
    igraph_time_goto(&graph, i);
    igraph_degree(&graph, &deg, igraph_vss_all(), IGRAPH_IN,
		  IGRAPH_NO_LOOPS);
    printf("  %i: ", i); igraph_vector_print(&deg);
  }

  igraph_vector_destroy(&deg);
  igraph_destroy(&graph);

  return 0;
}
