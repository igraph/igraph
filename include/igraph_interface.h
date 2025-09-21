/*
   igraph library.
   Copyright (C) 2009-2025  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef IGRAPH_INTERFACE_H
#define IGRAPH_INTERFACE_H

#include "igraph_decls.h"
#include "igraph_attributes.h"
#include "igraph_datatype.h"
#include "igraph_error.h"
#include "igraph_iterators.h"
#include "igraph_types.h"

IGRAPH_BEGIN_C_DECLS

/* -------------------------------------------------- */
/* Interface                                          */
/* -------------------------------------------------- */

IGRAPH_EXPORT igraph_error_t igraph_empty(igraph_t *graph, igraph_int_t n, igraph_bool_t directed);
IGRAPH_EXPORT igraph_error_t igraph_empty_attrs(
    igraph_t *graph, igraph_int_t n, igraph_bool_t directed,
    const igraph_attribute_record_list_t* attr
);
IGRAPH_EXPORT void igraph_destroy(igraph_t *graph);
IGRAPH_EXPORT igraph_error_t igraph_copy(igraph_t *to, const igraph_t *from);
IGRAPH_EXPORT igraph_error_t igraph_add_edges(
    igraph_t *graph, const igraph_vector_int_t *edges,
    const igraph_attribute_record_list_t* attr
);
IGRAPH_EXPORT igraph_error_t igraph_add_vertices(
    igraph_t *graph, igraph_int_t nv,
    const igraph_attribute_record_list_t* attr
);
IGRAPH_EXPORT igraph_error_t igraph_delete_edges(igraph_t *graph, igraph_es_t edges);
IGRAPH_EXPORT igraph_error_t igraph_delete_vertices(igraph_t *graph, igraph_vs_t vertices);
IGRAPH_EXPORT igraph_error_t igraph_delete_vertices_map(
    igraph_t *graph, igraph_vs_t vertices, igraph_vector_int_t *map,
    igraph_vector_int_t *invmap
);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_int_t igraph_vcount(const igraph_t *graph);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_int_t igraph_ecount(const igraph_t *graph);
IGRAPH_EXPORT igraph_error_t igraph_neighbors(
    const igraph_t *graph, igraph_vector_int_t *neis, igraph_int_t vid,
    igraph_neimode_t mode, igraph_loops_t loops, igraph_bool_t multiple
);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t igraph_is_directed(const igraph_t *graph);
IGRAPH_EXPORT igraph_error_t igraph_degree_1(
    const igraph_t *graph, igraph_int_t *deg, igraph_int_t vid,
    igraph_neimode_t mode, igraph_loops_t loops
);
IGRAPH_EXPORT igraph_error_t igraph_degree(
    const igraph_t *graph, igraph_vector_int_t *res, igraph_vs_t vids,
    igraph_neimode_t mode, igraph_loops_t loops
);
IGRAPH_EXPORT igraph_error_t igraph_edge(const igraph_t *graph, igraph_int_t eid,
                              igraph_int_t *from, igraph_int_t *to);
IGRAPH_EXPORT igraph_error_t igraph_edges(const igraph_t *graph, igraph_es_t eids,
                               igraph_vector_int_t *edges, igraph_bool_t bycol);
IGRAPH_EXPORT igraph_error_t igraph_get_eid(const igraph_t *graph, igraph_int_t *eid,
                                 igraph_int_t from, igraph_int_t to,
                                 igraph_bool_t directed, igraph_bool_t error);
IGRAPH_EXPORT igraph_error_t igraph_get_eids(const igraph_t *graph, igraph_vector_int_t *eids,
                                  const igraph_vector_int_t *pairs,
                                  igraph_bool_t directed, igraph_bool_t error);
IGRAPH_EXPORT igraph_error_t igraph_get_all_eids_between(const igraph_t *graph, igraph_vector_int_t *eids,
                                  igraph_int_t source, igraph_int_t target, igraph_bool_t directed);
IGRAPH_EXPORT igraph_error_t igraph_incident(
    const igraph_t *graph, igraph_vector_int_t *eids, igraph_int_t pnode,
    igraph_neimode_t mode, igraph_loops_t loops
);
IGRAPH_EXPORT igraph_error_t igraph_is_same_graph(const igraph_t *graph1, const igraph_t *graph2, igraph_bool_t *res);

IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t igraph_i_property_cache_get_bool(const igraph_t *graph, igraph_cached_property_t prop);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t igraph_i_property_cache_has(const igraph_t *graph, igraph_cached_property_t prop);
IGRAPH_EXPORT void igraph_i_property_cache_set_bool(const igraph_t *graph, igraph_cached_property_t prop, igraph_bool_t value);
IGRAPH_EXPORT void igraph_i_property_cache_set_bool_checked(const igraph_t *graph, igraph_cached_property_t prop, igraph_bool_t value);
IGRAPH_EXPORT void igraph_i_property_cache_invalidate(const igraph_t *graph, igraph_cached_property_t prop);
IGRAPH_EXPORT void igraph_i_property_cache_invalidate_all(const igraph_t *graph);

#define IGRAPH_RETURN_IF_CACHED_BOOL(graphptr, prop, resptr) \
    do { \
        if (igraph_i_property_cache_has((graphptr), (prop))) { \
            *(resptr) = igraph_i_property_cache_get_bool((graphptr), (prop)); \
            return IGRAPH_SUCCESS; \
        } \
    } while (0)

/**
 * \define IGRAPH_FROM
 * \brief The source vertex of an edge.
 *
 * Faster than \ref igraph_edge(), but no error checking is done: \p eid is assumed to be valid.
 *
 * \param graph The graph.
 * \param eid   The edge ID.
 * \return The source vertex of the edge.
 * \sa \ref igraph_edge() if error checking is desired.
 */
#define IGRAPH_FROM(graph,eid) ((igraph_int_t)(VECTOR((graph)->from)[(igraph_int_t)(eid)]))

/**
 * \define IGRAPH_TO
 * \brief The target vertex of an edge.
 *
 * Faster than \ref igraph_edge(), but no error checking is done: \p eid is assumed to be valid.
 *
 * \param graph The graph object.
 * \param eid   The edge ID.
 * \return The target vertex of the edge.
 * \sa \ref igraph_edge() if error checking is desired.
 */
#define IGRAPH_TO(graph,eid)   ((igraph_int_t)(VECTOR((graph)->to)  [(igraph_int_t)(eid)]))

/**
 * \define IGRAPH_OTHER
 * \brief The other endpoint of an edge.
 *
 * Typically used with undirected edges when one endpoint of the edge is known,
 * and the other endpoint is needed. No error checking is done:
 * \p eid and \p vid are assumed to be valid.
 *
 * \param graph The graph object.
 * \param eid   The edge ID.
 * \param vid   The vertex ID of one endpoint of an edge.
 * \return The other endpoint of the edge.
 * \sa \ref IGRAPH_TO() and \ref IGRAPH_FROM() to get the source and target
 *     of directed edges.
 */
#define IGRAPH_OTHER(graph,eid,vid) \
    ((igraph_int_t)(IGRAPH_TO(graph,(eid))==(vid) ? IGRAPH_FROM((graph),(eid)) : IGRAPH_TO((graph),(eid))))

IGRAPH_DEPRECATED IGRAPH_EXPORT igraph_error_t igraph_delete_vertices_idx(
    igraph_t *graph, igraph_vs_t vertices, igraph_vector_int_t *idx,
    igraph_vector_int_t *invidx
);

IGRAPH_END_C_DECLS

#endif
