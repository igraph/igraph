/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2006  Gabor Csardi <csardi@rmki.kfki.hu>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include "igraph.h"
#include "config.h"
#include <math.h>

/**
 * \function igraph_laplacian
 * \brief Returns the Laplacian matrix of a graph
 * 
 * </para><para>
 * The graph Laplacian matrix is similar to an adjacency matrix but
 * contains -1's instead of 1's and the vertex degrees are included in
 * the diagonal. So the result for edge i--j is -1 if i!=j and is equal
 * to the degree of vertex i if i==j. igraph_laplacian will work on a
 * directed graph; in this case, the diagonal will contain the out-degrees.
 * Loop edges will be ignored.
 * 
 * </para><para>
 * The normalized version of the Laplacian matrix has 1 in the diagonal and 
 * -1/sqrt(d[i]d[j]) if there is an edge from i to j.
 * 
 * </para><para>
 * The first version of this function was written by Vincent Matossian.
 * \param graph Pointer to the graph to convert.
 * \param res Pointer to an initialized matrix object, it will be
 *        resized if needed.
 * \param normalized Whether to create a normalized Laplacian matrix.
 * \return Error code.
 *
 * Time complexity: O(|V||V|),
 * |V| is the 
 * number of vertices in the graph.
 */

int igraph_laplacian(const igraph_t *graph, igraph_matrix_t *res,
                     igraph_bool_t normalized) {
  
  igraph_eit_t edgeit;
  long int no_of_nodes=igraph_vcount(graph);
  igraph_bool_t directed=igraph_is_directed(graph);
  long int from, to;
  igraph_integer_t ffrom, fto;
  igraph_vector_t degree;  
  int i;

  IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, no_of_nodes));
  igraph_matrix_null(res);
  IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(0), &edgeit));
  IGRAPH_FINALLY(igraph_eit_destroy, &edgeit);
  
  IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
  
  IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(),
                             IGRAPH_OUT, IGRAPH_NO_LOOPS));
  
  if(directed){
    if (!normalized) {
      for(i=0;i<no_of_nodes;i++)
        MATRIX(*res, i, i) = VECTOR(degree)[i];
      
      while (!IGRAPH_EIT_END(edgeit)) {
        igraph_edge(graph, IGRAPH_EIT_GET(edgeit), &ffrom, &fto);
        from=ffrom;
        to=fto;
        if (from != to) {
          MATRIX(*res, from, to) -= 1;
        }
        IGRAPH_EIT_NEXT(edgeit);
      }
    } else {
      for (i=0;i<no_of_nodes;i++) {
        MATRIX(*res, i, i) = VECTOR(degree)[i]>0 ? 1 : 0;
        if (VECTOR(degree)[i] > 0)
          VECTOR(degree)[i] = 1.0 / VECTOR(degree)[i];
      }
      
      while (!IGRAPH_EIT_END(edgeit)) {
        igraph_edge(graph, IGRAPH_EIT_GET(edgeit), &ffrom, &fto);
        from=ffrom; to=fto;
        if (from != to) {
          MATRIX(*res, from, to) -= VECTOR(degree)[from]; 
        }
        IGRAPH_EIT_NEXT(edgeit);
      }
    }

  } else {

    if (!normalized) {
      for(i=0;i<no_of_nodes;i++) {
        MATRIX(*res, i, i) = VECTOR(degree)[i];
      }
      
      while (!IGRAPH_EIT_END(edgeit)) {
        igraph_edge(graph, IGRAPH_EIT_GET(edgeit), &ffrom, &fto);
        from=ffrom;
        to=fto; 

        if (from != to) {
          MATRIX(*res, to, from) -= 1;
          MATRIX(*res, from, to) -= 1;
        }

        IGRAPH_EIT_NEXT(edgeit);
      }
    } else {
      for (i=0;i<no_of_nodes;i++) {
        MATRIX(*res, i, i) = VECTOR(degree)[i]>0 ? 1: 0;
        VECTOR(degree)[i] = sqrt(VECTOR(degree)[i]);
      }
      
      while (!IGRAPH_EIT_END(edgeit)) {
        igraph_edge(graph, IGRAPH_EIT_GET(edgeit), &ffrom, &fto);
        from=ffrom; to=fto;
        if (from != to) {
          double diff = 1.0 / (VECTOR(degree)[from] * VECTOR(degree)[to]);
          MATRIX(*res, from, to) -= diff;
          MATRIX(*res, to, from) -= diff;
        }
        IGRAPH_EIT_NEXT(edgeit);
      }
    }

  }

  igraph_vector_destroy(&degree);
  igraph_eit_destroy(&edgeit);
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}
