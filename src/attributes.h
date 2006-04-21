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

/* -------------------------------------------------- */
/* Attributes                                         */
/* -------------------------------------------------- */

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
 *    graph. The number of vertices to add is supplied as an
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
 */

typedef struct igraph_attribute_table_t {
  int (*init)(igraph_t *graph);
  void (*destroy)(igraph_t *graph);
  int (*copy)(igraph_t *to, const igraph_t *from);
  int (*add_vertices)(igraph_t *graph, long int nv);
  void (*delete_vertices)(igraph_t *graph, const igraph_vector_t *eidx,
			  const igraph_vector_t *vidx);
  int (*add_edges)(igraph_t *graph, long int ne);
  void (*delete_edges)(igraph_t *graph, const igraph_vector_t *idx);
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
        if (from->attr) { \
          IGRAPH_CHECK(igraph_i_ret=igraph_i_attribute_copy(to, from)); \
        } \
        if (igraph_i_ret != 0) { \
          IGRAPH_ERROR("", igraph_i_ret); \
        } \
   } while(0)        

int igraph_i_attribute_init(igraph_t *graph);
void igraph_i_attribute_destroy(igraph_t *graph);
int igraph_i_attribute_copy(igraph_t *to, const igraph_t *from);
int igraph_i_attribute_add_vertices(igraph_t *graph, long int nv);
void igraph_i_attribute_delete_vertices(igraph_t *graph, 
					const igraph_vector_t *eidx,
					const igraph_vector_t *vidx);
int igraph_i_attribute_add_edges(igraph_t *graph, long int ne);
void igraph_i_attribute_delete_edges(igraph_t *graph, 
				     const igraph_vector_t *idx);

#endif
