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

/**
 * \ingroup structural
 * \function igraph_are_connected
 * \brief Decides whether two vertices are connected 
 *
 * \param graph The graph object.
 * \param v1 The first vertex.
 * \param v2 The second vertex.
 * \return Boolean, \c TRUE if there is an edge from
 *         \p v1 to
 *         \p v2. Returns
 *         \c FALSE if there   
 *         is no such edge or at least one of the vertex ids is
 *         invalid (ie. too big or negative).
 * 
 * The function is of course symmetric for undirected graphs.
 *
 * Time complexity: O(d),
 * d is the
 * out-degree of \p v1.
 */
igraph_bool_t igraph_are_connected(const igraph_t *graph, 
			    igraph_integer_t v1, igraph_integer_t v2) {

  igraph_vs_t it;
  igraph_bool_t res=0;
  long int nov=igraph_vcount(graph);

  if (v1 < 0 || v2 < 0 || v1 > nov-1 || v2 > nov-1) {
    return 0;
/*     IGRAPH_ERROR("are connected", IGRAPH_EINVVID); */
  }

  IGRAPH_CHECK(igraph_vs_adj(graph, &it, v1, IGRAPH_OUT));
  
  while (!res && !igraph_vs_end(graph, &it)) {
    res= (igraph_vs_get(graph, &it) == v2);
    igraph_vs_next(graph, &it);
  }
  
  return res;
}
