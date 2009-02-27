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

#include "igraph.h"
#include "types.h"
#include "config.h"

/**
 * \ingroup structural
 * \function igraph_are_connected
 * \brief Decides whether two vertices are connected 
 *
 * \param graph The graph object.
 * \param v1 The first vertex.
 * \param v2 The second vertex.
 * \param res Boolean, \c TRUE if there is an edge from
 *         \p v1 to \p v2, \c FALSE otherwise.
 * \return Error code.
 * 
 * The function is of course symmetric for undirected graphs.
 *
 * </para><para>
 * Time complexity: O( min(log(d1), log(d2)) ),
 * d1 is the (out-)degree of \p v1 and d2 is the (in-)degree of \p v2.
 */
int igraph_are_connected(const igraph_t *graph, 
			 igraph_integer_t v1, igraph_integer_t v2,
			 igraph_bool_t *res) {

  long int nov=igraph_vcount(graph);
  igraph_integer_t eid=-1;

  if (v1 < 0 || v2 < 0 || v1 > nov-1 || v2 > nov-1) {
    IGRAPH_ERROR("are connected", IGRAPH_EINVVID);
  }

  igraph_get_eid2(graph, &eid, v1, v2, /*directed=*/1);
  *res = (eid >= 0);

  return IGRAPH_SUCCESS;
}
