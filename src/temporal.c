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

#include "igraph_temporal.h"
#include "igraph_interface.h"

int igraph_time_next(igraph_t *graph) {
  if (graph->now != IGRAPH_END) { graph->now += 1; }
  return 0;
}

int igraph_time_prev(igraph_t *graph) {
  if (graph->now > IGRAPH_BEGINNING) { graph->now -= 1; }
  return 0;
}

int igraph_time_goto(igraph_t *graph, igraph_time_t at) {
  graph->now = at;
  return 0;
}

int igraph_time_reset(igraph_t *graph) {
  graph->now = IGRAPH_BEGINNING;
  return 0;
}

int igraph_create_temporal(igraph_t *graph,
                           const igraph_vector_t *edges,
                           igraph_integer_t n, igraph_bool_t directed,
                           const igraph_vector_time_t *e_active,
                           const igraph_vector_time_t *e_inactive,
                           const igraph_vector_time_t *v_active,
                           const igraph_vector_time_t *v_inactive) {

  igraph_real_t rmin, rmax;
  igraph_integer_t min, max;
  igraph_integer_t no_verts=n;
  
  if (igraph_vector_size(edges) > 0) {
    igraph_vector_minmax(edges, &rmin, &rmax);
    min = (int) rmin;
    max = (int) rmax + 1;
  } else {
    min = max = 0;
  }
  
  if (igraph_vector_size(edges) % 2 != 0) {
    IGRAPH_ERROR("Invalid (odd) edges vector length",
                 IGRAPH_EINVEVECTOR);
  }

  if (min < 0) {
    IGRAPH_ERROR("Invalid (negative) vertex id", IGRAPH_EINVVID);
  }
  
  if (max > no_verts) { no_verts = max; }
  
  IGRAPH_CHECK(igraph_empty(graph, 0, directed));
  IGRAPH_FINALLY(igraph_destroy, graph);
  
  IGRAPH_CHECK(igraph_add_vertices_at(graph, no_verts, v_active,
                                      v_inactive, /*attr=*/ 0));
  IGRAPH_CHECK(igraph_add_edges_at(graph, edges, e_active,
                                   e_inactive, /*attr=*/ 0));

  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}
