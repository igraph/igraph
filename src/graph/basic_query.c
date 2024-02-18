/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "igraph_datatype.h"
#include "igraph_types.h"
#include "igraph_interface.h"
#include "igraph_structural.h"

/**
 * \ingroup structural
 * \function igraph_are_adjacent
 * \brief Decides whether two vertices are adjacent.
 *
 * Decides whether there are any edges that have \p v1 and \p v2
 * as endpoints. This function is of course symmetric for undirected
 * graphs.
 *
 * \param graph The graph object.
 * \param v1 The first vertex.
 * \param v2 The second vertex.
 * \param res Boolean, \c true if there is an edge from
 *         \p v1 to \p v2, \c false otherwise.
 * \return The error code \c IGRAPH_EINVVID is returned if an invalid
 *         vertex ID is given.
 *
 * Time complexity: O( min(log(d1), log(d2)) ),
 * d1 is the (out-)degree of \p v1 and d2 is the (in-)degree of \p v2.
 */
igraph_error_t igraph_are_adjacent(const igraph_t *graph,
                         igraph_integer_t v1, igraph_integer_t v2,
                         igraph_bool_t *res) {

    igraph_integer_t nov = igraph_vcount(graph);
    igraph_integer_t eid = -1;

    if (v1 < 0 || v2 < 0 || v1 > nov - 1 || v2 > nov - 1) {
        IGRAPH_ERROR("Invalid vertex ID when checking if two vertices are connected.", IGRAPH_EINVVID);
    }

    igraph_get_eid(graph, &eid, v1, v2, IGRAPH_DIRECTED, /*error=*/ false);
    *res = (eid >= 0);

    return IGRAPH_SUCCESS;
}


/**
 * \ingroup structural
 * \function igraph_are_connected
 * \brief Decides whether two vertices are adjacent (deprecated alias).
 *
 * \deprecated-by igraph_are_adjacent 0.10.10
 *
 * Decides whether there are any edges that have \p v1 and \p v2
 * as endpoints. This function is of course symmetric for undirected
 * graphs.
 *
 * \param graph The graph object.
 * \param v1 The first vertex.
 * \param v2 The second vertex.
 * \param res Boolean, \c true if there is an edge from
 *         \p v1 to \p v2, \c false otherwise.
 * \return The error code \c IGRAPH_EINVVID is returned if an invalid
 *         vertex ID is given.
 *
 * Time complexity: O( min(log(d1), log(d2)) ),
 * d1 is the (out-)degree of \p v1 and d2 is the (in-)degree of \p v2.
 */
igraph_error_t igraph_are_connected(const igraph_t *graph,
                                   igraph_integer_t v1, igraph_integer_t v2,
                                   igraph_bool_t *res) {
    return igraph_are_adjacent(graph, v1, v2, res);
}
