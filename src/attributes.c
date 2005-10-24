/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2005  Gabor Csardi <csardi@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
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
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/

#include "igraph.h"
#include "memory.h"

int igraph_add_graph_attribute(igraph_t *graph, const char *name) {
  igraph_attribute_list_add(&graph->gal, name);
  return 0;
}

int igraph_remove_graph_attribute(igraph_t *graph, const char *name) {
  igraph_attribute_list_remove(&graph->gal, name);  
  return 0;
}

int igraph_get_graph_attribute(igraph_t *graph, const char *name, 
			       real_t *value) {
  igraph_attribute_list_get(&graph->gal, name, 0, value);
  return 0;
}

int igraph_set_graph_attribute(igraph_t *graph, const char *name, 
			       real_t value) {
  igraph_attribute_list_set(&graph->gal, name, 0, value);
  return 0;
}

int igraph_list_graph_attributes(igraph_t *graph, igraph_strarray_t *l) {
  igraph_attribute_list_list(&graph->gal, l);
  return 0;
}

int igraph_add_vertex_attribute(igraph_t *graph, const char *name) {
  igraph_attribute_list_add(&graph->val, name);
  return 0;
}

int igraph_remove_vertex_attribute(igraph_t *graph, const char *name) {
  igraph_attribute_list_remove(&graph->val, name);
  return 0;
}

int igraph_get_vertex_attribute(igraph_t *graph, const char *name, 
				long int v, real_t *value) {
  igraph_attribute_list_get(&graph->val, name, v, value);
  return 0;  
}

int igraph_set_vertex_attribute(igraph_t *graph, const char *name, 
				long int v, real_t value) {
  igraph_attribute_list_set(&graph->val, name, v, value);
  return 0;
}

int igraph_get_vertex_attributes(igraph_t *graph, const char *name, 
				 vector_t *v, vector_t *value) {
  igraph_attribute_list_gets(&graph->val, name, v, value);
  return 0;  
}

int igraph_set_vertex_attributes(igraph_t *graph, const char *name, 
				 vector_t *v, vector_t *value) {
  igraph_attribute_list_sets(&graph->val, name, v, value);
  return 0;
}

int igraph_list_vertex_attributes(igraph_t *graph, igraph_strarray_t *l) {
  igraph_attribute_list_list(&graph->val, l);
  return 0;
}
