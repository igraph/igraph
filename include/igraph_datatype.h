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

#ifndef IGRAPH_DATATYPE_H
#define IGRAPH_DATATYPE_H

#include "igraph_decls.h"
#include "igraph_types.h"
#include "igraph_vector.h"

__BEGIN_DECLS

struct igraph_i_property_cache_t;
typedef struct igraph_i_property_cache_t igraph_i_property_cache_t;

typedef enum {
    /* Stores whether the graph has at least one self-loop. */
    IGRAPH_PROP_HAS_LOOP = 0,

    /* Stores whether the graph has at least one multi-edge, taking into account
     * edge directions in directed graphs. In other words, this property should
     * be false for a directed graph with edges (a, b) and (b, a), and true
     * for a directed graph with edges (a, b) and (a, b) again. */
    IGRAPH_PROP_HAS_MULTI,

    /* Stores whether the graph has at least one reciprocal edge pair. Ignored
     * in undirected graphs. This property should be true for a directed graph
     * with edges (a, b) and (b, a), and false for a directed graph with
     * edges (a, b) and (a, b) again. Self-loops (a, a) are not considered
     * reciprocal. */
    IGRAPH_PROP_HAS_MUTUAL,

    /* Stores whether the graph is weakly connected. */
    IGRAPH_PROP_IS_WEAKLY_CONNECTED,

    /* Stores whether the graph is strongly connected. Ignored in undirected graphs. */
    IGRAPH_PROP_IS_STRONGLY_CONNECTED,

    /* Stores whether the graph is a directed acyclic graph. Not used for
     * undirected graphs. */
    IGRAPH_PROP_IS_DAG,

    /* Stores whether the graph is a forest, i.e. an undirected or directed
     * graph that is cycle-free even if we ignore edge directions. */
    IGRAPH_PROP_IS_FOREST,

    /* Dummy value used to count enum values */
    IGRAPH_PROP_I_SIZE
} igraph_cached_property_t;

/**
 * \ingroup internal
 * \struct igraph_t
 * \brief The internal data structure for storing graphs.
 *
 * It is simple and efficient. It has the following members:
 * - <b>n</b> The number of vertices, redundant.
 * - <b>directed</b> Whether the graph is directed.
 * - <b>from</b> The first column of the edge list.
 * - <b>to</b> The second column of the edge list.
 * - <b>oi</b> The index of the edge list by the first column. Thus
 *   the first edge according to this order goes from
 *   \c from[oi[0]] to \c to[oi[0]]. The length of
 *   this vector is the same as the number of edges in the graph.
 * - <b>ii</b> The index of the edge list by the second column.
 *   The length of this vector is the same as the number of edges.
 * - <b>os</b> Contains pointers to the edgelist (\c from
 *   and \c to for every vertex. The first edge \em from
 *   vertex \c v is edge no. \c from[oi[os[v]]] if
 *   \c os[v]<os[v+1]. If \c os[v]==os[v+1] then
 *   there are no edges \em from node \c v. Its length is
 *   the number of vertices plus one, the last element is always the
 *   same as the number of edges and is contained only to ease the
 *   queries.
 * - <b>is</b> This is basically the same as <b>os</b>, but this time
 *   for the incoming edges.
 *
 * For undirected graphs, the same edge list is stored, i.e. an
 * undirected edge is stored only once. Currently, undirected edges
 * are canonicalized so that the index of the 'from' vertex is not greater
 * than the index of the 'to' vertex. Thus, if v1 <= v2, only the edge (v1, v2)
 * needs to be searched for, not (v2, v1), to determine if v1 and v2 are connected.
 * However, this fact is NOT guaranteed by the documented public API,
 * and should not be relied upon by the implementation of any functions,
 * except those belonging to the minimal API in type_indexededgelist.c.
 *
 * The storage requirements for a graph with \c |V| vertices
 * and \c |E| edges is \c O(|E|+|V|).
 */
typedef struct igraph_s {
    igraph_integer_t n;
    igraph_bool_t directed;
    igraph_vector_int_t from;
    igraph_vector_int_t to;
    igraph_vector_int_t oi;
    igraph_vector_int_t ii;
    igraph_vector_int_t os;
    igraph_vector_int_t is;
    void *attr;
    igraph_i_property_cache_t *cache;
} igraph_t;

IGRAPH_EXPORT void igraph_invalidate_cache(const igraph_t* graph);

__END_DECLS

#endif
