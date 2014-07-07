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
  return igraph_data_type_temp_time_next
    ((igraph_data_type_temp_t *) graph);
}

int igraph_time_prev(igraph_t *graph) {
  return igraph_data_type_temp_time_prev
    ((igraph_data_type_temp_t *) graph);
  return 0;
}

int igraph_time_goto(igraph_t *graph, igraph_time_t at) {
  return igraph_data_type_temp_time_goto
    ((igraph_data_type_temp_t *) graph, at);
  return 0;
}

int igraph_time_reset(igraph_t *graph) {
  return igraph_data_type_temp_time_reset
    ((igraph_data_type_temp_t *) graph);
  return 0;
}

igraph_time_t igraph_now(const igraph_t *graph) {
  return igraph_data_type_temp_now((igraph_data_type_temp_t *) graph);
}

int igraph_time_last(const igraph_t *graph, igraph_time_t *vertex,
                     igraph_time_t *edge,
                     igraph_time_t *vertex_or_edge) {

  igraph_time_t v, e;

  if (vertex || vertex_or_edge) {
    if (igraph_vector_int_is_null(&graph->vb)) {
      v = IGRAPH_END;
    } else {
      v = igraph_vector_int_size(&graph->vb) - 2;
    }
    if (vertex) { *vertex = v; }
  }
  if (edge || vertex_or_edge) {
    if (igraph_vector_int_is_null(&graph->eb)) {
      e = IGRAPH_END;
    } else {
      e = igraph_vector_int_size(&graph->eb) - 2;
    }
    if (edge) { *edge = e; }
  }
  if (vertex_or_edge) {
    *vertex_or_edge = v > e ? v : e;
  }
  return 0;
}

int igraph_vertices_range(const igraph_t *graph, igraph_vs_t vs,
			  igraph_vector_time_t *active,
			  igraph_vector_time_t *inactive) {

  /* Cannot do vcount here... */
  igraph_integer_t vs_nodes;
  igraph_vs_size(graph, &vs, &vs_nodes);

  if (active) {
    IGRAPH_CHECK(igraph_vector_time_resize(active, vs_nodes));
    if (igraph_vector_int_is_null(&graph->vb)) {
      igraph_vector_time_null(active);
    } else {
      if (igraph_vs_is_all(&vs)) {
	int i, j = 0, l = igraph_vector_int_size(&graph->vb) - 1;
	for (i = 0; i < l; i++) {
	  int f = VECTOR(graph->vb)[i];
	  int t = VECTOR(graph->vb)[i + 1];
	  for ( ; f < t; f++) {
	    VECTOR(*active)[j++] = i;
	  }
	}
      } else {
	igraph_vit_t vit;
	int i = 0;
	igraph_vit_create(graph, vs, &vit);
	for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit);
	     IGRAPH_VIT_NEXT(vit)) {
	  igraph_integer_t node = IGRAPH_VIT_GET(vit);
	  long int pos;
	  igraph_bool_t found =
	    igraph_vector_int_binsearch(&graph->vb, node, &pos);
	  VECTOR(*active)[i++] = found ? pos : pos - 1;
	}
      }
    }
  }

  if (inactive) {
    IGRAPH_CHECK(igraph_vector_time_resize(inactive, vs_nodes));
    if (igraph_vector_time_is_null(&graph->vd)) {
      igraph_vector_time_fill(inactive, IGRAPH_END);
    } else {
      if (igraph_vs_is_all(&vs)) {
	igraph_vector_time_update(inactive, &graph->vd);
      } else {
	igraph_vit_t vit;
	int i = 0;
	igraph_vit_create(graph, vs, &vit);
	for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit);
	     IGRAPH_VIT_NEXT(vit)) {
	  igraph_integer_t node = IGRAPH_VIT_GET(vit);
	  VECTOR(*inactive)[i++] = VECTOR(graph->vd)[node];
	}
      }
    }
  }

  return 0;
}

int igraph_edges_range(const igraph_t *graph, igraph_es_t es,
		       igraph_vector_time_t *active,
		       igraph_vector_time_t *inactive) {

  /* Cannot do ecount here... */
  igraph_integer_t es_edges;
  igraph_es_size(graph, &es, &es_edges);

  if (active) {
    IGRAPH_CHECK(igraph_vector_time_resize(active, es_edges));
    if (igraph_vector_int_is_null(&graph->eb)) {
      igraph_vector_time_null(active);
    } else {
      if (igraph_es_is_all(&es)) {
	int i, j = 0, l = igraph_vector_int_size(&graph->eb) - 1;
	for (i = 0; i < l; i++) {
	  int f = VECTOR(graph->eb)[i];
	  int t = VECTOR(graph->eb)[i + 1];
	  for ( ; f < t; f++) {
	    VECTOR(*active)[j++] = i;
	  }
	}
      } else {
	igraph_eit_t eit;
	int i = 0;
	igraph_eit_create(graph, es, &eit);
	for (IGRAPH_EIT_RESET(eit); !IGRAPH_EIT_END(eit);
	     IGRAPH_EIT_NEXT(eit)) {
	  igraph_integer_t edge = IGRAPH_EIT_GET(eit);
	  long int pos;
	  igraph_bool_t found =
	    igraph_vector_int_binsearch(&graph->eb, edge, &pos);
	  VECTOR(*active)[i++] = found ? pos : pos - 1;
	}
      }
    }
  }

  if (inactive) {
    IGRAPH_CHECK(igraph_vector_time_resize(inactive, es_edges));
    if (igraph_vector_time_is_null(&graph->ed)) {
      igraph_vector_time_fill(inactive, IGRAPH_END);
    } else {
      if (igraph_es_is_all(&es)) {
	igraph_vector_time_update(inactive, &graph->ed);
      } else {
	igraph_eit_t eit;
	int i = 0;
	igraph_eit_create(graph, es, &eit);
	for (IGRAPH_EIT_RESET(eit); !IGRAPH_EIT_END(eit);
	     IGRAPH_EIT_NEXT(eit)) {
	  igraph_integer_t edge=IGRAPH_EIT_GET(eit);
	  VECTOR(*inactive)[i++] = VECTOR(graph->ed)[edge];
	}
      }
    }
  }

  return 0;
}

int igraph_time_slice(const igraph_t *graph, igraph_t *result,
                      igraph_time_t from, igraph_time_t to) {

  igraph_integer_t v_from, v_to, e_from, e_to;
  igraph_vector_t del_edges, del_vertices;
  int no_del_edges, no_del_vertices;
  igraph_time_t last_v, last_e, last;
  int no_nodes = graph->n, no_edges = igraph_vector_size(&graph->from);
  int i, j;

  igraph_time_last(graph, &last_v, &last_e, &last);

  if (from == IGRAPH_END) { from = last; }
  if (to   == IGRAPH_END) { to   = last; }

  v_from = from == IGRAPH_BEGINNING ? 0 : VECTOR(graph->vb)[from - 1];
  v_to   = to >= last_v ? no_nodes - 1 : VECTOR(graph->vb)[to + 1] - 1;

  e_from = from == IGRAPH_BEGINNING ? 0 : VECTOR(graph->eb)[from - 1];
  e_to   = to >= last_e ? no_edges -1 : VECTOR(graph->eb)[to + 1] - 1;

  IGRAPH_CHECK(igraph_copy(result, graph));
  IGRAPH_FINALLY(igraph_destroy, result);

  no_del_edges = no_edges - (e_to - e_from + 1);
  if (no_del_edges > 0) {
    IGRAPH_VECTOR_INIT_FINALLY(&del_edges, no_del_edges);
    for (i = 0, j = 0; i < e_from  ; i++) { VECTOR(del_edges)[j++] = i; }
    for (i = e_to + 1; i < no_edges; i++) { VECTOR(del_edges)[j++] = i; }

    IGRAPH_CHECK(igraph_delete_edges(result,
				     igraph_ess_vector(&del_edges)));
    igraph_vector_destroy(&del_edges);
    IGRAPH_FINALLY_CLEAN(1);
  }

  no_del_vertices = no_nodes - (v_to - v_from + 1);
  if (no_del_vertices > 0) {
    IGRAPH_VECTOR_INIT_FINALLY(&del_vertices, no_del_vertices);
    for (i = 0, j = 0; i < v_from  ; i++) { VECTOR(del_vertices)[j++] = i; }
    for (i = v_to + 1; i < no_nodes; i++) { VECTOR(del_vertices)[j++] = i; }

    IGRAPH_CHECK(igraph_delete_vertices(result,
					igraph_vss_vector(&del_vertices)));
    igraph_vector_destroy(&del_vertices);
    IGRAPH_FINALLY_CLEAN(1);
  }

  IGRAPH_FINALLY_CLEAN(1); 	/* result */

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
