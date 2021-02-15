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

#ifndef IGRAPH_VISITOR_H
#define IGRAPH_VISITOR_H

#include "igraph_decls.h"
#include "igraph_constants.h"
#include "igraph_types.h"
#include "igraph_datatype.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Visitor-like functions                             */
/* -------------------------------------------------- */

/**
 * \typedef igraph_bfshandler_t
 * Callback type for BFS function
 *
 * \ref igraph_bfs() is able to call a callback function, whenever a
 * new vertex is found, while doing the breadth-first search. This
 * callback function must be of type \c igraph_bfshandler_t. It has
 * the following arguments:
 * \param graph The graph that that algorithm is working on. Of course
 *   this must not be modified.
 * \param vid The id of the vertex just found by the breadth-first
 *   search.
 * \param pred The id of the previous vertex visited. It is -1 if
 *   there is no previous vertex, because the current vertex is the root
 *   is a search tree.
 * \param succ The id of the next vertex that will be visited. It is
 *   -1 if there is no next vertex, because the current vertex is the
 *   last one in a search tree.
 * \param rank The rank of the current vertex, it starts with zero.
 * \param dist The distance (number of hops) of the current vertex
 *   from the root of the current search tree.
 * \param extra The extra argument that was passed to \ref
 *   igraph_bfs().
 * \return A logical value, if TRUE (=non-zero), that is interpreted
 *    as a request to stop the BFS and return to the caller. If a BFS
 *    is terminated like this, then all elements of the result vectors
 *    that were not yet calculated at the point of the termination
 *    contain NaN.
 *
 * \sa \ref igraph_bfs()
 */

typedef igraph_bool_t igraph_bfshandler_t(const igraph_t *graph,
        igraph_integer_t vid,
        igraph_integer_t pred,
        igraph_integer_t succ,
        igraph_integer_t rank,
        igraph_integer_t dist,
        void *extra);

IGRAPH_EXPORT int igraph_bfs(const igraph_t *graph,
                             igraph_integer_t root, const igraph_vector_t *roots,
                             igraph_neimode_t mode, igraph_bool_t unreachable,
                             const igraph_vector_t *restricted,
                             igraph_vector_t *order, igraph_vector_t *rank,
                             igraph_vector_t *father,
                             igraph_vector_t *pred, igraph_vector_t *succ,
                             igraph_vector_t *dist, igraph_bfshandler_t *callback,
                             void *extra);

IGRAPH_EXPORT int igraph_bfs_simple(igraph_t *graph, igraph_integer_t vid, igraph_neimode_t mode,
                                    igraph_vector_t *vids, igraph_vector_t *layers,
                                    igraph_vector_t *parents);

/**
 * \function igraph_dfshandler_t
 * Callback type for the DFS function
 *
 * \ref igraph_dfs() is able to call a callback function, whenever a
 * new vertex is discovered, and/or whenever a subtree is
 * completed. These callbacks must be of type \c
 * igraph_dfshandler_t. They have the following arguments:
 * \param graph The graph that that algorithm is working on. Of course
 *   this must not be modified.
 * \param vid The id of the vertex just found by the depth-first
 *   search.
 * \param dist The distance (number of hops) of the current vertex
 *   from the root of the current search tree.
 * \param extra The extra argument that was passed to \ref
 *   igraph_dfs().
 * \return A logical value, if TRUE (=non-zero), that is interpreted
 *    as a request to stop the DFS and return to the caller. If a DFS
 *    is terminated like this, then all elements of the result vectors
 *    that were not yet calculated at the point of the termination
 *    contain NaN.
 *
 * \sa \ref igraph_dfs()
 */

typedef igraph_bool_t igraph_dfshandler_t(const igraph_t *graph,
        igraph_integer_t vid,
        igraph_integer_t dist,
        void *extra);

IGRAPH_EXPORT int igraph_dfs(const igraph_t *graph, igraph_integer_t root,
                             igraph_neimode_t mode, igraph_bool_t unreachable,
                             igraph_vector_t *order,
                             igraph_vector_t *order_out, igraph_vector_t *father,
                             igraph_vector_t *dist, igraph_dfshandler_t *in_callback,
                             igraph_dfshandler_t *out_callback,
                             void *extra);

__END_DECLS

#endif
