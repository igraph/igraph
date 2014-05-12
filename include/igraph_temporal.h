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

#ifndef IGRAPH_TEMPORAL_H
#define IGRAPH_TEMPORAL_H

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
#include "igraph_vector.h"
#include "igraph_datatype.h"

__BEGIN_DECLS

#ifdef IGRAPH_DATA_TYPE_TEMPORAL_EDGE_LIST

/* ---------------------------------------------------------------- */

/**
 * \function igraph_empty_at
 */
 
int igraph_empty_at(igraph_t *graph, igraph_integer_t n,
                igraph_bool_t directed,
                const igraph_vector_time_t *v_active,
                const igraph_vector_time_t *v_inactive,
                void *attrs);

igraph_integer_t igraph_vcount_at(const igraph_t *graph,
                igraph_time_t at);
igraph_integer_t igraph_ecount_at(const igraph_t *graph,
                igraph_time_t at);

int igraph_neighbors_at(const igraph_t *graph, igraph_vector_t *neis,
                igraph_integer_t vid, igraph_neimode_t mode,
                igraph_time_t at);

int igraph_incident_at(const igraph_t *graph, igraph_vector_t *eids,
                igraph_integer_t vid, igraph_neimode_t mode,
                igraph_time_t at);
                      
int igraph_degree_at(const igraph_t *graph, igraph_vector_t *res, 
		const igraph_vs_t vids, igraph_neimode_t mode, 
		igraph_bool_t loops, igraph_time_t at);
int igraph_get_eid_at(const igraph_t *graph, igraph_integer_t *eid,
		igraph_integer_t from, igraph_integer_t to,
		igraph_bool_t directed, igraph_bool_t error,
                igraph_time_t at);
int igraph_get_eids_at(const igraph_t *graph, igraph_vector_t *eids,
                const igraph_vector_t *pairs,
		const igraph_vector_t *path,
		igraph_bool_t directed, igraph_bool_t error,
                igraph_time_t at);
int igraph_get_eids_multi_at(const igraph_t *graph,
                igraph_vector_t *eids, const igraph_vector_t *pairs, 
	        const igraph_vector_t *path, igraph_bool_t directed,
                igraph_bool_t error, igraph_time_t at);

int igraph_add_edges_at(igraph_t *graph, const igraph_vector_t *edges,
                const igraph_vector_int_t *e_active,
                const igraph_vector_int_t *e_inactive, void *attr);
int igraph_add_vertices_at(igraph_t *graph, igraph_integer_t nv, 
		const igraph_vector_int_t *v_active,
                const igraph_vector_int_t *v_inacive, void *attr);

/* ---------------------------------------------------------------- */

/**
 * \function igraph_create_temporal
 *
 * Time steps are numbered from zero.
 * \param graph An uninitialized graph object.
 * \param edges The edges to add, the first two elements are the first
 *        edge, etc.
 * \param n The number of vertices in the graph, if smaller or equal
 *        than the highest vertex id in the \p edges vector it
 *        will be increased automatically. So it is safe to give 0
 *        here.
 * \param directed Boolean, whether to create a directed graph or
 *        not. If yes, then the first edge points from the first
 *        vertex id in \p edges to the second, etc.
 * \param e_active Vector the specifying the times steps in which the 
 *        edges become active in the network. If this is a null pointer,
 *        then all edges are active since the beginning of time (time
 *        step zero).
 * \param e_inactive Vector giving the time steps in which the edges
 *        become inactive in the network. The vector must be
 *        non-decreasing, i.e. it is expected that edges are ordered
 *        according to their activation. If edges live forever, then
 *        this can be a null pointer. If only a subset of the edges
 *        live forever, then set their values to -1. It is OK to include
 *        values that are larger than <code>nt</code>.
 * \param v_active Integer vector, giving the time step when each vertex
 *        becomes active. This vector must be non-decreasing, in other
 *        words vertices must be ordered accorcing to their activation
 *        time. If this is a null pointer, then all vertices
 *        are active since the beginning of time.
 * \param v_inactive Integer vertor, giving the time steps when each
 *        vertex becomes inactive. If all vertices live forever, then
 *        supply a null pointer here. If only a subset of vertices live
 *        forever, then supply -1 for them.
 * \return Error code.
 *
 */

int igraph_create_temporal(igraph_t *graph,
                           const igraph_vector_t *edges,
                           igraph_integer_t n, igraph_bool_t directed,
                           const igraph_vector_int_t e_active,
                           const igraph_vector_int_t e_inactive,
                           const igraph_vector_int_t v_active,
                           const igraph_vector_int_t v_inactive);

#endif

__END_DECLS

#endif
