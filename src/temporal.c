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

/* -------------------------------------------------- */
/* Interface, temporal functions                      */
/* -------------------------------------------------- */

int igraph_time_next(igraph_t *graph) {
  graph->now += 1;
  return 0;
}

int igraph_time_prev(igraph_t *graph) {
  if (graph->now != 0) { graph->now -= 1; }
  return 0;
}

int igraph_time_goto(igraph_t *graph, igraph_time_t at) {
  if (at != IGRAPH_END && at < IGRAPH_BEGINNING) {
    IGRAPH_ERROR("Invalid time step, cannot set time cursor",
                  IGRAPH_EINVAL);
  }
  graph->now = at;
  return 0;
}

int igraph_time_reset(igraph_t *graph) {
  graph->now = IGRAPH_BEGINNING;
}

int igraph_add_edges_at(igraph_t *graph, const igraph_vector_t *edges,
                const igraph_vector_int_t *e_active,
                const igraph_vector_int_t *e_inactive, void *attr) {
  /* TODO */
}

int igraph_add_vertices_at(igraph_t *graph, igraph_integer_t nv, 
		const igraph_vector_int_t *v_active,
                const igraph_vector_int_t *v_inacive, void *attr) {
  /* TODO */
}

int igraph_create_temporal(igraph_t *graph,
                           const igraph_vector_t *edges,
                           igraph_integer_t n, igraph_bool_t directed,
                           const igraph_vector_time_t e_active,
                           const igraph_vector_time_t e_inactive,
                           const igraph_vector_time_t v_active,
                           const igraph_vector_time_t v_inactive) {
  /* TODO */
}
