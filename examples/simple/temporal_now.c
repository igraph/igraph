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

  IGRAPH_VECTOR_CONSTANT(edges, 0,1,1,2,2,3,3,4,4,5,5,0);

  igraph_create_temporal(&graph, &edges, 0, IGRAPH_DIRECTED, 0, 0, 0, 0);

  if (igraph_now(&graph) != IGRAPH_END) { return 1; }

  igraph_time_next(&graph);
  if (igraph_now(&graph) != IGRAPH_END) { return 2; }

  igraph_time_reset(&graph);
  if (igraph_now(&graph) != IGRAPH_BEGINNING) { return 3; }

  igraph_time_prev(&graph);
  if (igraph_now(&graph) != IGRAPH_BEGINNING) { return 4; }

  igraph_time_next(&graph);
  if (igraph_now(&graph) != 1) { return 5; }

  igraph_time_prev(&graph);
  if (igraph_now(&graph) != IGRAPH_BEGINNING) { return 6; }

  igraph_destroy(&graph);

  return 0;
}
