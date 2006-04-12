/* -*- mode: C -*-  */
/* 
   IGraph R package.
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
#include "memory.h"

int igraph_cocitation_real(const igraph_t *graph, igraph_matrix_t *res, 
			   const igraph_vs_t *vids, 
			   igraph_neimode_t mode);

/**
 * \ingroup structural
 * \function igraph_cocitation
 * \brief Cocitation coupling.
 * 
 * Two vertices are cocited if there is another vertex citing both of
 * them. \ref igraph_cocitation() simply counts how many types two vertices are
 * cocited.
 * The cocitation score for each given vertex and all other vertices
 * in the graph will be calculated.
 * \param graph The graph object to analyze.
 * \param res Pointer to a matrix, the result of the calculation will
 *        be stored here. The number of its rows is the same as the
 *        number of vertex ids in \p vids, the number of
 *        columns is the number of vertices in the graph.
 * \param vids The vertex ids of the vertices for which the
 *        calculation will be done.
 * \return Error code:
 *         \c IGRAPH_EINVVID: invalid vertex id.
 * 
 * Time complexity: O(|V|d^2), |V| is
 * the number of vertices in the graph,
 * d is the (maximum) degree of 
 * the vertices in the graph.
 *
 * \sa \ref igraph_bibcoupling()
 */

int igraph_cocitation(const igraph_t *graph, igraph_matrix_t *res, 
		      const igraph_vs_t *vids) {
  return igraph_cocitation_real(graph, res, vids, IGRAPH_OUT);
}

/**
 * \ingroup structural
 * \function igraph_bibcoupling
 * \brief Bibliographic coupling.
 * 
 * The bibliographic coupling of two vertices is the number
 * of other vertices they both cite, \ref igraph_bibcoupling() calculates
 * this.
 * The bibliographic coupling  score for each given vertex and all
 * other vertices in the graph will be calculated.
 * \param graph The graph object to analyze.
 * \param res Pointer to a matrix, the result of the calculation will
 *        be stored here. The number of its rows is the same as the
 *        number of vertex ids in \p vids, the number of
 *        columns is the number of vertices in the graph.
 * \param vids The vertex ids of the vertices for which the
 *        calculation will be done.
 * \return Error code:
 *         \c IGRAPH_EINVVID: invalid vertex id.
 * 
 * Time complexity: O(|V|d^2),
 * |V| is the number of vertices in
 * the graph, d is the (maximum)
 * degree of the vertices in the graph.
 *
 * \sa \ref igraph_cocitation()
 */

int igraph_bibcoupling(const igraph_t *graph, igraph_matrix_t *res, 
		       const igraph_vs_t *vids) {
  return igraph_cocitation_real(graph, res, vids, IGRAPH_IN);
}

int igraph_cocitation_real(const igraph_t *graph, igraph_matrix_t *res, 
			   const igraph_vs_t *vids,
			   igraph_neimode_t mode) {

  long int no_of_nodes=igraph_vcount(graph);
  long int from, i, j;
  igraph_bool_t *calc;
  igraph_matrix_t tmpres=IGRAPH_MATRIX_NULL;
  igraph_vector_t neis=IGRAPH_VECTOR_NULL;
  igraph_vs_t myvids;
  const igraph_vector_t *myvidsv;
  
  IGRAPH_CHECK(igraph_vs_vectorview_it(graph, vids, &myvids));
  IGRAPH_FINALLY(igraph_vs_destroy, &myvids);
  myvidsv=igraph_vs_vector_getvector(graph, &myvids);

  if (!igraph_vector_isininterval(myvidsv, 0, no_of_nodes-1)) {
    IGRAPH_ERROR("", IGRAPH_EINVVID);
  }
  
  calc=Calloc(no_of_nodes, igraph_bool_t);
  if (calc==0) {
    IGRAPH_ERROR("cannot calculate cocitation/bibcoupling", IGRAPH_ENOMEM);
  }  
  IGRAPH_FINALLY(free, calc); 	/* TODO: hack */

  for (i=0; i<igraph_vector_size(myvidsv); i++) {
    calc[ (long int) VECTOR(*myvidsv)[i] ] = 1;
  }
  
  IGRAPH_MATRIX_INIT_FINALLY(&tmpres, no_of_nodes, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_CHECK(igraph_matrix_resize(res, igraph_vector_size(myvidsv), no_of_nodes));

  /* The result */
  
  for (from=0; from<no_of_nodes; from++) {
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, from, mode));
    for (i=0; i < igraph_vector_size(&neis)-1; i++) {
      if (calc[ (long int)VECTOR(neis)[i] ]) {
	for (j=i+1; j<igraph_vector_size(&neis); j++) {
	  MATRIX(tmpres, (long int)VECTOR(neis)[i], 
		 (long int)VECTOR(neis)[j]) += 1;
	  MATRIX(tmpres, (long int)VECTOR(neis)[j], 
		 (long int)VECTOR(neis)[i]) += 1;
	}
      }
    }
  }

  /* Copy result */
  for (i=0; i<igraph_vector_size(myvidsv); i++) {
    for (j=0; j<no_of_nodes; j++) {
      MATRIX(*res, i, j) = MATRIX(tmpres, (long int) VECTOR(*myvidsv)[i], j);
    }
  }  
  
  /* Clean up */
  igraph_matrix_destroy(&tmpres);
  igraph_vector_destroy(&neis);
  Free(calc);
  igraph_vs_destroy(&myvids);
  IGRAPH_FINALLY_CLEAN(4);

  if (vids->shorthand) { igraph_vs_destroy((igraph_vs_t*)vids); }  
  return 0;
}
