/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#ifndef IGRAPH_INTERFACE_H
#define IGRAPH_INTERFACE_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

#include "igraph_types.h"
#include "igraph_datatype.h"
#include "igraph_iterators.h"

__BEGIN_DECLS

#ifdef IGRAPH_DATA_TYPE_INDEXED_EDGE_LIST
#  define igraph_empty               igraph_empty_ie
#  define igraph_empty_attrs         igraph_empty_attrs_ie
#  define igraph_destroy             igraph_destroy_ie
#  define igraph_copy                igraph_copy_ie
#  define igraph_add_edges           igraph_add_edges_ie
#  define igraph_add_vertices        igraph_add_vertices_ie
#  define igraph_delete_edges        igraph_delete_edges_ie
#  define igraph_delete_vertices     igraph_delete_vertices_ie
#  define igraph_delete_vertices_idx igraph_delete_vertices_idx_ie
#  define igraph_vcount              igraph_vcount_ie
#  define igraph_ecount              igraph_ecount_ie
#  define igraph_neighbors           igraph_neighbors_ie
#  define igraph_is_directed         igraph_is_directed_ie
#  define igraph_degree              igraph_degree_ie
#  define igraph_edge                igraph_edge_ie
#  define igraph_edges               igraph_edges_ie
#  define igraph_get_eid             igraph_get_eid_ie
#  define igraph_get_eids            igraph_get_eids_ie
#  define igraph_get_eids_multi      igraph_get_eids_multi_ie
#  define igraph_adjacent            igraph_adjacenct_ie
#  define igraph_incident            igraph_incident_ie
#  define IGRAPH_FROM                IGRAPH_FROM_IE
#  define IGRAPH_TO                  IGRAPH_TO_IE
#  define IGRAPH_OTHER               IGRAPH_OTHER_IE
#endif

#ifdef IGRAPH_DATA_TYPE_TEMPORAL_EDGE_LIST
#  define igraph_empty               igraph_empty_temp
#  define igraph_empty_attrs         igraph_empty_attrs_temp
#  define igraph_destroy             igraph_destroy_temp
#  define igraph_copy                igraph_copy_temp
#  define igraph_add_edges           igraph_add_edges_temp
#  define igraph_add_vertices        igraph_add_vertices_temp
#  define igraph_delete_edges        igraph_delete_edges_temp
#  define igraph_delete_vertices     igraph_delete_vertices_temp
#  define igraph_delete_vertices_idx igraph_delete_vertices_idx_temp
#  define igraph_vcount              igraph_vcount_temp
#  define igraph_ecount              igraph_ecount_temp
#  define igraph_neighbors           igraph_neighbors_temp
#  define igraph_is_directed         igraph_is_directed_temp
#  define igraph_degree              igraph_degree_temp
#  define igraph_edge                igraph_edge_temp
#  define igraph_edges               igraph_edges_temp
#  define igraph_get_eid             igraph_get_eid_temp
#  define igraph_get_eids            igraph_get_eids_temp
#  define igraph_get_eids_multi      igraph_get_eids_multi_temp
#  define igraph_adjacent            igraph_adjacenct_temp
#  define igraph_incident            igraph_incident_temp
#  define IGRAPH_FROM                IGRAPH_FROM_TEMP
#  define IGRAPH_TO                  IGRAPH_TO_TEMP
#  define IGRAPH_OTHER               IGRAPH_OTHER_TEMP
#endif

/* -------------------------------------------------- */
/* Interface                                          */
/* -------------------------------------------------- */

int igraph_empty_ie(igraph_data_type_ie_t *graph, igraph_integer_t n,
                    igraph_bool_t directed);
int igraph_empty_attrs_ie(igraph_data_type_ie_t *graph, igraph_integer_t n,
                          igraph_bool_t directed, void *attr);
int igraph_destroy_ie(igraph_data_type_ie_t *graph);
int igraph_copy_ie(igraph_data_type_ie_t *to, const igraph_data_type_ie_t *from);
int igraph_add_edges_ie(igraph_data_type_ie_t *graph, const igraph_vector_t *edges, 
		        void *attr);
int igraph_add_vertices_ie(igraph_data_type_ie_t *graph, igraph_integer_t nv, 
			   void *attr);
int igraph_delete_edges_ie(igraph_data_type_ie_t *graph, igraph_es_t edges);
int igraph_delete_vertices_ie(igraph_data_type_ie_t *graph,
                              const igraph_vs_t vertices);
int igraph_delete_vertices_idx_ie(igraph_data_type_ie_t *graph,
                                  const igraph_vs_t vertices, 
			          igraph_vector_t *idx, 
			          igraph_vector_t *invidx);
igraph_integer_t igraph_vcount_ie(const igraph_data_type_ie_t *graph);
igraph_integer_t igraph_ecount_ie(const igraph_data_type_ie_t *graph);
int igraph_neighbors_ie(const igraph_data_type_ie_t *graph, igraph_vector_t *neis,
                        igraph_integer_t vid, igraph_neimode_t mode); 
igraph_bool_t igraph_is_directed_ie(const igraph_data_type_ie_t *graph);
int igraph_degree_ie(const igraph_data_type_ie_t *graph, igraph_vector_t *res, 
		     const igraph_vs_t vids, igraph_neimode_t mode, 
		     igraph_bool_t loops);
int igraph_edge_ie(const igraph_data_type_ie_t *graph, igraph_integer_t eid, 
	   	igraph_integer_t *from, igraph_integer_t *to);		
int igraph_edges_ie(const igraph_data_type_ie_t *graph, igraph_es_t eids,
		    igraph_vector_t *edges);
int igraph_get_eid_ie(const igraph_data_type_ie_t *graph, igraph_integer_t *eid,
		    igraph_integer_t from, igraph_integer_t to,
		    igraph_bool_t directed, igraph_bool_t error);
int igraph_get_eids_ie(const igraph_data_type_ie_t *graph, igraph_vector_t *eids,
                       const igraph_vector_t *pairs,
		       const igraph_vector_t *path,
	               igraph_bool_t directed, igraph_bool_t error);
int igraph_get_eids_multi_ie(const igraph_data_type_ie_t *graph,
                             igraph_vector_t *eids,
			     const igraph_vector_t *pairs, 
			     const igraph_vector_t *path,
			     igraph_bool_t directed,
                             igraph_bool_t error);
/* deprecated */
int igraph_adjacent_ie(const igraph_data_type_ie_t *graph, igraph_vector_t *eids,
                       igraph_integer_t vid, igraph_neimode_t mode);
int igraph_incident_ie(const igraph_data_type_ie_t *graph, igraph_vector_t *eids,
                       igraph_integer_t vid, igraph_neimode_t mode);

#define IGRAPH_FROM_IE(g,e) \
        ((igraph_integer_t)(VECTOR((g)->from)[(long int)(e)]))
#define IGRAPH_TO_IE(g,e)   \
        ((igraph_integer_t)(VECTOR((g)->to)  [(long int)(e)]))
#define IGRAPH_OTHER_IE(g,e,v) \
        ((igraph_integer_t)(IGRAPH_TO(g,(e))==(v) ? \
          IGRAPH_FROM((g),(e)) : IGRAPH_TO((g),(e))))

/* ----------------------------------------------------------------- */

int igraph_empty_temp(igraph_data_type_temp_t *graph, igraph_integer_t n,
                    igraph_bool_t directed);
int igraph_empty_attrs_temp(igraph_data_type_temp_t *graph, igraph_integer_t n,
                          igraph_bool_t directed, void *attr);
int igraph_destroy_temp(igraph_data_type_temp_t *graph);
int igraph_copy_temp(igraph_data_type_temp_t *to, const igraph_data_type_temp_t *from);
int igraph_add_edges_temp(igraph_data_type_temp_t *graph, const igraph_vector_t *edges, 
		        void *attr);
int igraph_add_vertices_temp(igraph_data_type_temp_t *graph, igraph_integer_t nv, 
			   void *attr);
int igraph_delete_edges_temp(igraph_data_type_temp_t *graph, igraph_es_t edges);
int igraph_delete_vertices_temp(igraph_data_type_temp_t *graph,
                              const igraph_vs_t vertices);
int igraph_delete_vertices_idx_temp(igraph_data_type_temp_t *graph,
                                  const igraph_vs_t vertices, 
			          igraph_vector_t *idx, 
			          igraph_vector_t *invidx);
igraph_integer_t igraph_vcount_temp(const igraph_data_type_temp_t *graph);
igraph_integer_t igraph_ecount_temp(const igraph_data_type_temp_t *graph);
int igraph_neighbors_temp(const igraph_data_type_temp_t *graph, igraph_vector_t *neis,
                        igraph_integer_t vid, igraph_neimode_t mode); 
igraph_bool_t igraph_is_directed_temp(const igraph_data_type_temp_t *graph);
int igraph_degree_temp(const igraph_data_type_temp_t *graph, igraph_vector_t *res, 
		     const igraph_vs_t vids, igraph_neimode_t mode, 
		     igraph_bool_t loops);
int igraph_edge_temp(const igraph_data_type_temp_t *graph, igraph_integer_t eid, 
	   	igraph_integer_t *from, igraph_integer_t *to);		
int igraph_edges_temp(const igraph_data_type_temp_t *graph, igraph_es_t eids,
		    igraph_vector_t *edges);
int igraph_get_eid_temp(const igraph_data_type_temp_t *graph, igraph_integer_t *eid,
		    igraph_integer_t from, igraph_integer_t to,
		    igraph_bool_t directed, igraph_bool_t error);
int igraph_get_eids_temp(const igraph_data_type_temp_t *graph, igraph_vector_t *eids,
                       const igraph_vector_t *pairs,
		       const igraph_vector_t *path,
	               igraph_bool_t directed, igraph_bool_t error);
int igraph_get_eids_multi_temp(const igraph_data_type_temp_t *graph,
                             igraph_vector_t *eids,
			     const igraph_vector_t *pairs, 
			     const igraph_vector_t *path,
			     igraph_bool_t directed,
                             igraph_bool_t error);
/* deprecated */
int igraph_adjacent_temp(const igraph_data_type_temp_t *graph, igraph_vector_t *eids,
                       igraph_integer_t vid, igraph_neimode_t mode);
int igraph_incident_temp(const igraph_data_type_temp_t *graph, igraph_vector_t *eids,
                       igraph_integer_t vid, igraph_neimode_t mode);

#define IGRAPH_FROM_TEMP(g,e) \
        ((igraph_integer_t)(VECTOR((g)->from)[(long int)(e)]))
#define IGRAPH_TO_TEMP(g,e)   \
        ((igraph_integer_t)(VECTOR((g)->to)  [(long int)(e)]))
#define IGRAPH_OTHER_TEMP(g,e,v) \
        ((igraph_integer_t)(IGRAPH_TO(g,(e))==(v) ? \
          IGRAPH_FROM((g),(e)) : IGRAPH_TO((g),(e))))

/**
 * \function igraph_add_edges_at
 * Add edges to a temporal graph
 *
 * TODO
 */

int igraph_add_edges_at(igraph_data_type_temp_t *graph,
		const igraph_vector_t *edges,
                const igraph_vector_time_t *e_active,
                const igraph_vector_time_t *e_inactive, void *attr);

/**
 * \function igraph_add_vertices_at
 * Add vertices to a temporal graph
 */
                
int igraph_add_vertices_at(igraph_data_type_temp_t *graph,
	        igraph_integer_t nv, const igraph_vector_time_t *v_active,
                const igraph_vector_time_t *v_inacive, void *attr);


__END_DECLS

#endif
