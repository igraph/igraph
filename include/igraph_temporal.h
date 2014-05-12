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

/**
 * \section temporal_graphs Temporal graphs
 * 
 * Temporal graphs' structrue changes over time: edges are deleted
 * and/or added, vertices are deleted and/or added over time.
 *
 * igraph stores the complete time evolution of temporal graphs,
 * and makes it possible to run queries on graph, and any given
 * point of time.
 *
 * Time steps are denoted by non-negative integers. 
 * We assume that the first time step (the beginning of time) is
 * time step 0 (\ref IGRAPH_BEGINNING), and the last time step is
 * the end of time (\ref IGRAPH_END).
 * 
 * We say that a vertex or edge was activated in a given time step \c t
 * if it was not present in the previous time step, and it is present
 * in \c t. We say that a vertex or edge was inactivated in a given time
 * step \c t, if it was present in the previous time step, and it is not
 * present in \c t.
 *
 * Vertices and edges cannot be activated again, once they were
 * inactivated. (But of course it is possible to activate _another_
 * edge with the same pair of incident vertices.)
 *
 * 
 * TODO: now
 */
 
/* ---------------------------------------------------------------- */

/**
 * \function igraph_time_next
 * Slide the time cursor to the next time step
 * 
 * TODO
 */

int igraph_time_next(igraph_t *graph);

/**
 * \function igraph_time_prev
 * Slide the time cursor to the previous time step
 *
 * TODO
 */

int igraph_time_prev(igraph_t *graph);

/**
 * \function igraph_time_goto
 * Set the time cursor to the given time step
 *
 * TODO
 */

int igraph_time_goto(igraph_t *graph, igraph_time_t at);

/**
 * \function igraph_time_reset
 * Set the time cursor to the beginning of time
 * 
 * TODO
 */

int igraph_time_reset(igraph_t *graph);

/**
 * \function igraph_add_edges_at
 * Add edges to a temporal graph
 *
 * TODO
 */

int igraph_add_edges_at(igraph_t *graph, const igraph_vector_t *edges,
                const igraph_vector_int_t *e_active,
                const igraph_vector_int_t *e_inactive, void *attr);

/**
 * \function igraph_add_vertices_at
 * Add vertices to a temporal graph
 */
                
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
 *        live forever, then set their values to \c IGRAPH_END. 
 * \param v_active Vector, giving the time step when each vertex
 *        becomes active. This vector must be non-decreasing, in other
 *        words vertices must be ordered accorcing to their activation
 *        time. If this is a null pointer, then all vertices
 *        are active since the beginning of time.
 * \param v_inactive Vector, giving the time steps when each
 *        vertex becomes inactive. If all vertices live forever, then
 *        supply a null pointer here. If only a subset of vertices live
 *        forever, then supply \c IGRAPH_END for them.
 * \return Error code.
 *
 */

int igraph_create_temporal(igraph_t *graph,
                           const igraph_vector_t *edges,
                           igraph_integer_t n, igraph_bool_t directed,
                           const igraph_vector_time_t e_active,
                           const igraph_vector_time_t e_inactive,
                           const igraph_vector_time_t v_active,
                           const igraph_vector_time_t v_inactive);

#endif

__END_DECLS

#endif
