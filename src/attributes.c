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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include "attributes.h"
#include "memory.h"
#include "igraph.h"

/**
 * \section about_attributes
 * 
 * <para>Attributes are numbers or strings (or basically any data
 * structure) associated with the vertices or edges of a graph, or
 * with the graph itself. Eg. you may associate symbolic names with
 * vertices or numeric weights with the edges of a graph. </para>
 * 
 * <para>Attribute handling has been largely changed in \a igraph
 * 0.2. From now on it is possible to attach an attribute handling
 * interface to \a igraph. This is simply a table of functions, of
 * type \ref igraph_attribute_table_t. This functions are invoked to
 * notify the attribute handling code about the structural changes in
 * a graph. See the documentation of this type for details.</para>
 *
 * <para>By default there is no attribute interface attached to \a igraph,
 * to attach one, call \ref igraph_i_set_attribute_table with your new
 * table. </para>
 *
 */

int igraph_i_attribute_init(igraph_t *graph, void *attr) {
  graph->attr=0;
  if (igraph_i_attribute_table) {
    return igraph_i_attribute_table->init(graph, attr);
  } else {
    return 0;
  }
}
  
void igraph_i_attribute_destroy(igraph_t *graph) {
  if (igraph_i_attribute_table) {
    igraph_i_attribute_table->destroy(graph);
  }
}
  
int igraph_i_attribute_copy(igraph_t *to, const igraph_t *from) {
  if (igraph_i_attribute_table) {
    return igraph_i_attribute_table->copy(to, from);
  } else {
    return 0;
    }
  }
  
int igraph_i_attribute_add_vertices(igraph_t *graph, long int nv, void *attr) {
  if (igraph_i_attribute_table) {
    return igraph_i_attribute_table->add_vertices(graph, nv, attr);
  } else {
    return 0;
  }
}

void igraph_i_attribute_delete_vertices(igraph_t *graph, 
					const igraph_vector_t *eidx,
					const igraph_vector_t *vidx) {
  if (igraph_i_attribute_table) {
    igraph_i_attribute_table->delete_vertices(graph, eidx, vidx);
  }
}
  
  
int igraph_i_attribute_add_edges(igraph_t *graph, 
				 const igraph_vector_t *edges, void *attr) {
  if (igraph_i_attribute_table) {
    return igraph_i_attribute_table->add_edges(graph, edges, attr);
  } else { 
    return 0;
  }
}
  
void igraph_i_attribute_delete_edges(igraph_t *graph, 
				     const igraph_vector_t *idx) {
  if (igraph_i_attribute_table) {
    igraph_i_attribute_table->delete_edges(graph, idx);
  }
}

int igraph_i_attribute_permute_edges(igraph_t *graph, 
				      const igraph_vector_t *idx) {
  if (igraph_i_attribute_table) {
    return igraph_i_attribute_table->permute_edges(graph, idx);
  } else {
    return 0;
  }
}

int igraph_i_attribute_get_info(const igraph_t *graph,
				igraph_strvector_t *gnames, 
				igraph_vector_t *gtypes,
				igraph_strvector_t *vnames,
				igraph_vector_t *vtypes,
				igraph_strvector_t *enames,
				igraph_vector_t *etypes) {
  if (igraph_i_attribute_table) {
    return igraph_i_attribute_table->get_info(graph, gnames, gtypes, 
	  			              vnames, vtypes, 
				              enames, etypes);
  } else {
    return 0;
  }
}

igraph_bool_t igraph_i_attribute_has_attr(const igraph_t *graph, 
					  igraph_attribute_elemtype_t type,
					  const char *name) {
  if (igraph_i_attribute_table) {
    return igraph_i_attribute_table->has_attr(graph, type, name);
  } else {
    return 0;
  }
}

int igraph_i_attribute_gettype(const igraph_t *graph,
			       igraph_attribute_type_t *type,
			       igraph_attribute_elemtype_t elemtype,
			       const char *name) {
  if (igraph_i_attribute_table) {
    return igraph_i_attribute_table->gettype(graph, type, elemtype, name);
  } else {
    return 0;
  }
  
}

int igraph_i_attribute_get_numeric_graph_attr(const igraph_t *graph,
					      const char *name,
					      igraph_vector_t *value) {
  if (igraph_i_attribute_table) {
    return igraph_i_attribute_table->get_numeric_graph_attr(graph, name, value);
  } else {
    return 0;
  }
}

int igraph_i_attribute_get_numeric_vertex_attr(const igraph_t *graph, 
					       const char *name,
					       igraph_vs_t vs,
					       igraph_vector_t *value) {
  if (igraph_i_attribute_table) {
    return igraph_i_attribute_table->get_numeric_vertex_attr(graph, name, vs, value);
  } else {
    return 0;
  }
}

int igraph_i_attribute_get_numeric_edge_attr(const igraph_t *graph,
					     const char *name,
					     igraph_es_t es,
					     igraph_vector_t *value) {
  if (igraph_i_attribute_table) {
    return igraph_i_attribute_table->get_numeric_edge_attr(graph, name, es, value);
  } else {
    return 0;
  }
}

int igraph_i_attribute_get_string_graph_attr(const igraph_t *graph,
					     const char *name,
					     igraph_strvector_t *value) {
  if (igraph_i_attribute_table) {
    return igraph_i_attribute_table->get_string_graph_attr(graph, name, value);
  } else {
    return 0;
  }
}

int igraph_i_attribute_get_string_vertex_attr(const igraph_t *graph, 
					      const char *name,
					      igraph_vs_t vs,
					      igraph_strvector_t *value) {
  if (igraph_i_attribute_table) {
    return igraph_i_attribute_table->get_string_vertex_attr(graph, name, vs, value);
  } else {
    return 0;
  }
}

int igraph_i_attribute_get_string_edge_attr(const igraph_t *graph,
					    const char *name,
					    igraph_es_t es,
					    igraph_strvector_t *value) {
  if (igraph_i_attribute_table) {
    return igraph_i_attribute_table->get_string_edge_attr(graph, name, es, value);
  } else {
    return 0;
  }
}

igraph_attribute_table_t *igraph_i_attribute_table=0;

/**
 * \function igraph_i_set_attribute_table
 * This function attaches attribute handling code to the igraph library.
 * \param table Pointer to an \ref igraph_attribute_table_t object
 *    containing the functions for attribute manipulation. Supply \c
 *    NULL here if you don't want attributes.
 * \return Pointer to the old attribute handling table.
 * 
 * Time complexity: O(1).
 */
  
igraph_attribute_table_t *
igraph_i_set_attribute_table(igraph_attribute_table_t * table) {
  igraph_attribute_table_t *old=igraph_i_attribute_table;
  igraph_i_attribute_table=table;
  return old;
}
  
