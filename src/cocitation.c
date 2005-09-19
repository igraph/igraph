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
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/

#include "igraph.h"
#include "memory.h"

int igraph_cocitation_real(igraph_t *graph, matrix_t *res, vector_t *vids, 
			   igraph_neimode_t mode);

/**
 * \ingroup structural
 * \brief Cocitation coupling
 * 
 * Two vertices are cocited if there is another vertex citing both of
 * them. igraph_cocitation() siply counts how many types two vertices are
 * cocited. The bibliographic coupling of two vertices is the number
 * of other vertices they both cite, igraph_bibcoupling() calculates
 * this.
 * The cocitation score for each given vertex and all other vertices
 * in the graph will be calculated.
 * @param graph The graph object to analyze.
 * @param res Pointer to a matrix, the result of the calculation will
 *        be stored here. The number of its rows is the same as the
 *        number of vertex ids in <code>vids</code>, the number of
 *        columns is the number of vertices in the graph.
 * @param vids The vertex ids of the vertices for which the
 *        calculation will be done.
 * @return Error code.
 * 
 * Time complexity: <code>O(|V|d^2</code>, |V| is the number of
 * vertices in the graph, <code>d</code> is the (maximum) degree of
 * the vertices in the graph.
 *
 * \sa igraph_bibcoupling()
 */

int igraph_cocitation(igraph_t *graph, matrix_t *res, vector_t *vids) {
  return igraph_cocitation_real(graph, res, vids, IGRAPH_OUT);
}

/**
 * \ingroup structural
 * \brief Bibliographic coupling
 * 
 * Two vertices are cocited if there is another vertex citing both of
 * them. igraph_cocitation() siply counts how many types two vertices are
 * cocited. The bibliographic coupling of two vertices is the number
 * of other vertices they both cite, igraph_bibcoupling() calculates
 * this.
 * The bibliographic coupling  score for each given vertex and all
 * other vertices in the graph will be calculated.
 * @param graph The graph object to analyze.
 * @param res Pointer to a matrix, the result of the calculation will
 *        be stored here. The number of its rows is the same as the
 *        number of vertex ids in <code>vids</code>, the number of
 *        columns is the number of vertices in the graph.
 * @param vids The vertex ids of the vertices for which the
 *        calculation will be done.
 * @return Error code.
 * 
 * Time complexity: <code>O(|V|d^2</code>, |V| is the number of
 * vertices in the graph, <code>d</code> is the (maximum) degree of
 * the vertices in the graph.
 *
 * \sa igraph_cocitation()
 */

int igraph_bibcoupling(igraph_t *graph, matrix_t *res, vector_t *vids) {
  return igraph_cocitation_real(graph, res, vids, IGRAPH_IN);
}

int igraph_cocitation_real(igraph_t *graph, matrix_t *res, vector_t *vids, 
			   igraph_neimode_t mode) {

  long int no_of_nodes=igraph_vcount(graph);
  long int from, i, j;
  bool_t *calc;
  matrix_t tmpres;
  vector_t neis;
  
  calc=Calloc(no_of_nodes, bool_t);
  for (i=0; i<vector_size(vids); i++) {
    calc[ (long int) VECTOR(*vids)[i] ] = 1;
  }
  
  matrix_init(&tmpres, no_of_nodes, no_of_nodes);
  vector_init(&neis, 0);

  /* The result */
  
  for (from=0; from<no_of_nodes; from++) {
    igraph_neighbors(graph, &neis, from, mode);
    for (i=0; i < vector_size(&neis)-1; i++) {
      if (calc[ (long int)VECTOR(neis)[i] ]) {
	for (j=i+1; j<vector_size(&neis); j++) {
	  MATRIX(tmpres, (long int)VECTOR(neis)[i], 
		 (long int)VECTOR(neis)[j]) += 1;
	  MATRIX(tmpres, (long int)VECTOR(neis)[j], 
		 (long int)VECTOR(neis)[i]) += 1;
	}
      }
    }
  }

  /* Copy result */
  matrix_resize(res, vector_size(vids), no_of_nodes);
  for (i=0; i<vector_size(vids); i++) {
    for (j=0; j<no_of_nodes; j++) {
      MATRIX(*res, i, j) = MATRIX(tmpres, (long int) VECTOR(*vids)[i], j);
    }
  }  
  
  /* Clean up */
  matrix_destroy(&tmpres);
  vector_destroy(&neis);
  
  return 0;
}
