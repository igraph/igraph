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
 * \param res Boolean, \c TRUE if there is an edge from
 *         \p v1 to \p v2, \c FALSE otherwise.
 * \return Error code.
 * 
 * The function is of course symmetric for undirected graphs.
 *
 * </para><para>
 * Time complexity: O(d),
 * d is the
 * out-degree of \p v1.
 */
int igraph_are_connected(const igraph_t *graph, 
			 igraph_integer_t v1, igraph_integer_t v2,
			 igraph_bool_t *res) {
  igraph_vs_t vs;
  igraph_vit_t it;
  long int nov=igraph_vcount(graph);

  if (v1 < 0 || v2 < 0 || v1 > nov-1 || v2 > nov-1) {
    IGRAPH_ERROR("are connected", IGRAPH_EINVVID);
  }

  IGRAPH_CHECK(igraph_vs_adj(&vs, v1, IGRAPH_OUT));
  IGRAPH_CHECK(igraph_vit_create(graph, vs, &it));
  IGRAPH_FINALLY(igraph_vit_destroy, &it);
  
  *res=0;
  IGRAPH_VIT_RESET(it);
  while (!*res && !IGRAPH_VIT_END(it)) {
    *res= (IGRAPH_VIT_GET(it) == v2);
    IGRAPH_VIT_NEXT(it);
  }
  
  igraph_vit_destroy(&it);
  IGRAPH_FINALLY_CLEAN(1);
  
  return IGRAPH_SUCCESS;
}
