/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2010  Gabor Csardi <csardi.gabor@gmail.com>
   Rue de l'Industrie 5, Lausanne 1005, Switzerland
   
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include <igraph.h>
#include <stdio.h>

#define FAIL(msg, error) do { printf(msg "\n") ; return error; } while (0)

int main() {
  
  igraph_t graph;
  igraph_vector_long_t sep;
  igraph_bool_t result;

  igraph_vector_long_init(&sep, 0);

  /* Simple star graph, remove the center */
  igraph_star(&graph, 10, IGRAPH_STAR_UNDIRECTED, 0);
  igraph_vector_long_resize(&sep, 1);
  VECTOR(sep)[0] = 0;
  igraph_is_minimal_separator(&graph, &sep, &result);
  if (!result) FAIL("Center of star graph failed.", 1);

  /* Same graph, but another vertex */
  VECTOR(sep)[0] = 6;
  igraph_is_minimal_separator(&graph, &sep, &result);
  if (result) FAIL("Non-center of star graph failed.", 2);
  igraph_destroy(&graph);

  /* Karate club */
  igraph_famous(&graph, "zachary");
  igraph_vector_long_clear(&sep);
  igraph_vector_long_push_back(&sep, 32);
  igraph_vector_long_push_back(&sep, 33);
  igraph_is_minimal_separator(&graph, &sep, &result);
  if (!result) FAIL("Karate network (32,33) failed", 3);

  igraph_vector_long_resize(&sep, 5);
  VECTOR(sep)[0]=8;
  VECTOR(sep)[1]=9;
  VECTOR(sep)[2]=19;
  VECTOR(sep)[3]=30;
  VECTOR(sep)[4]=31;
  igraph_is_minimal_separator(&graph, &sep, &result);
  if (result) FAIL("Karate network (8,9,19,30,31) failed", 4);

  igraph_destroy(&graph);
  igraph_vector_long_destroy(&sep);

  return 0;
}

