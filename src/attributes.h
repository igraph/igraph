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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#ifndef REST_ATTRIBUTES_H
#define REST_ATTRIBUTES_H

#include "igraph.h"

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Attributes                                         */
/* -------------------------------------------------- */

typedef enum { IGRAPH_ATTRIBUTE_DEFAULT=0,
	       IGRAPH_ATTRIBUTE_NUMERIC=1, 
	       IGRAPH_ATTRIBUTE_STRING=2,
	       IGRAPH_ATTRIBUTE_R_OBJECT=3, 
	       IGRAPH_ATTRIBUTE_PY_OBJECT=4 } igraph_attribute_type_t;

typedef struct igraph_i_attribute_record_t {
  const char *name;
  igraph_attribute_type_t type;
  const void *value;
} igraph_i_attribute_record_t;

typedef enum { IGRAPH_ATTRIBUTE_GRAPH=0,
	       IGRAPH_ATTRIBUTE_VERTEX,
	       IGRAPH_ATTRIBUTE_EDGE } igraph_attribute_elemtype_t;

/**
 * \struct igraph_attribute_table_t
 * 
 * This type collects the functions defining an attribute handler.
 * It has the following members:
 * \member init This function is called whenever a new graph object is
 *    created, right after it is created but before any vertices or
 *    edges are added. It is supposed to set the \c attr member of the \c
 *    igraph_t object. It is expected to return an error code.
 * \member destroy This function is called whenever the graph object
 *    is destroyed, right before freeing the allocated memory. 
 * \member copy This function is called when copying a graph with \ref
 *    igraph_copy, after the structure of the graph has been already
 *    copied. It is expected to return an error code.
 * \member add_vertices Called when vertices are added to a
 *    graph, before adding the vertices themselves.
 *    The number of vertices to add is supplied as an
 *    argument. Expected to return an error code. 
 * \member delete_vertices Called when vertices are deleted from the
 *    graph. Two additional parameters are supplied, the first is a
 *    recoding vector for edge ids, the second is one for the vertex
 *    ids. The edge recoding vector gives for each edge its id in the
 *    new graph. It contains one number for each edge (in the original
 *    graph): zero means that the edge has been deleted, otherwise the
 *    new id plus one is included. The vertex recoding vector contains
 *    the same for vertices.
 * \member add_edges Called when new edges have been added. The number
 *    of new edges are supplied as well. It is expected to return an
 *    error code.
 * \member delete_edges Called when edges were deleted. The edge
 *    recoding vector is supplied, in the same form as for the \c
 *    delete_vertices function.
 * \member permute_edges Typically called when a new graph is created and 
 *    some of the new edges should carry the attributes of some of the
 *    old edges. The idx vector shows the mapping between the old edges and 
 *    the new ones. Its length is the same as the number of edges in the new 
 *    graph, and for each edge it gives the id of the old edge (the edge in
 *    the old graph).
 */

typedef struct igraph_attribute_table_t {
  int (*init)(igraph_t *graph);
  void (*destroy)(igraph_t *graph);
  int (*copy)(igraph_t *to, const igraph_t *from);
  int (*add_vertices)(igraph_t *graph, long int nv, igraph_vector_ptr_t *attr);
  void (*delete_vertices)(igraph_t *graph, const igraph_vector_t *eidx,
			  const igraph_vector_t *vidx);
  int (*add_edges)(igraph_t *graph, const igraph_vector_t *edges, 
		   igraph_vector_ptr_t *attr);
  void (*delete_edges)(igraph_t *graph, const igraph_vector_t *idx);
  int (*permute_edges)(igraph_t *graph, const igraph_vector_t *idx);

  int (*get_info)(const igraph_t *graph,
		  igraph_strvector_t *gnames, igraph_vector_t *gtypes,
		  igraph_strvector_t *vnames, igraph_vector_t *vtypes,
		  igraph_strvector_t *enames, igraph_vector_t *etypes);

  igraph_bool_t (*has_attr)(const igraph_t *graph, igraph_attribute_elemtype_t type,
			    const char *name);
  int (*gettype)(const igraph_t *graph, igraph_attribute_type_t *type,
		 igraph_attribute_elemtype_t elemtype, const char *name);

  int (*get_numeric_graph_attr)(const igraph_t *graph, const char *name,
				igraph_vector_t *value);
  int (*get_string_graph_attr)(const igraph_t *graph, const char *name,
			       igraph_strvector_t *value);
  int (*get_numeric_vertex_attr)(const igraph_t *graph, const char *name,
				 igraph_vs_t vs,
				 igraph_vector_t *value);
  int (*get_string_vertex_attr)(const igraph_t *graph, const char *name,
				igraph_vs_t vs,
				igraph_strvector_t *value);
  int (*get_numeric_edge_attr)(const igraph_t *graph, const char *name,
			       igraph_es_t es,
			       igraph_vector_t *value);
  int (*get_string_edge_attr)(const igraph_t *graph, const char *name,
			      igraph_es_t es,
			      igraph_strvector_t *value);
} igraph_attribute_table_t;

extern igraph_attribute_table_t *igraph_i_attribute_table;

igraph_attribute_table_t *
igraph_i_set_attribute_table(igraph_attribute_table_t * table);

#define IGRAPH_I_ATTRIBUTE_DESTROY(graph) \
        do {if ((graph)->attr) igraph_i_attribute_destroy(graph);} while(0)
#define IGRAPH_I_ATTRIBUTE_DELETE_VERTICES(graph, eidx, vidx) \
        do {if ((graph)->attr) igraph_i_attribute_delete_vertices((graph),(eidx),(vidx));} while(0)
#define IGRAPH_I_ATTRIBUTE_COPY(to,from) do { \
        int igraph_i_ret=0; \
        if ((from)->attr) { \
          IGRAPH_CHECK(igraph_i_ret=igraph_i_attribute_copy((to), (from))); \
        } \
        if (igraph_i_ret != 0) { \
          IGRAPH_ERROR("", igraph_i_ret); \
        } \
   } while(0)        

int igraph_i_attribute_init(igraph_t *graph);
void igraph_i_attribute_destroy(igraph_t *graph);
int igraph_i_attribute_copy(igraph_t *to, const igraph_t *from);
int igraph_i_attribute_add_vertices(igraph_t *graph, long int nv, void *attr);
void igraph_i_attribute_delete_vertices(igraph_t *graph, 
					const igraph_vector_t *eidx,
					const igraph_vector_t *vidx);
int igraph_i_attribute_add_edges(igraph_t *graph, 
				 const igraph_vector_t *edges, void *attr);
void igraph_i_attribute_delete_edges(igraph_t *graph, 
				     const igraph_vector_t *idx);
int igraph_i_attribute_permute_edges(igraph_t *graph, 
				     const igraph_vector_t *idx);

int igraph_i_attribute_get_info(const igraph_t *graph,
				igraph_strvector_t *gnames, 
				igraph_vector_t *gtypes,
				igraph_strvector_t *vnames,
				igraph_vector_t *vtypes,
				igraph_strvector_t *enames,
				igraph_vector_t *etypes);
igraph_bool_t igraph_i_attribute_has_attr(const igraph_t *graph, 
					  igraph_attribute_elemtype_t type,
					  const char *name);
int igraph_i_attribute_gettype(const igraph_t *graph,
			       igraph_attribute_type_t *type,
			       igraph_attribute_elemtype_t elemtype,
			       const char *name);

int igraph_i_attribute_get_numeric_graph_attr(const igraph_t *graph,
					      const char *name,
					      igraph_vector_t *value);
int igraph_i_attribute_get_numeric_vertex_attr(const igraph_t *graph, 
					       const char *name,
					       igraph_vs_t vs,
					       igraph_vector_t *value);
int igraph_i_attribute_get_numeric_edge_attr(const igraph_t *graph,
					     const char *name,
					     igraph_es_t es,
					     igraph_vector_t *value);
int igraph_i_attribute_get_string_graph_attr(const igraph_t *graph,
					     const char *name,
					     igraph_strvector_t *value);
int igraph_i_attribute_get_string_vertex_attr(const igraph_t *graph, 
					      const char *name,
					      igraph_vs_t vs,
					      igraph_strvector_t *value);
int igraph_i_attribute_get_string_edge_attr(const igraph_t *graph,
					    const char *name,
					    igraph_es_t es,
					    igraph_strvector_t *value);

__END_DECLS

#endif
