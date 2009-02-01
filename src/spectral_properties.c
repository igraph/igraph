/* -*- mode: C -*-  */
/* vim:set ts=8 sw=2 sts=2 et: */
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

int igraph_i_weighted_laplacian(const igraph_t *graph, igraph_matrix_t *res,
				igraph_bool_t normalized, 
				igraph_laplacian_direction_t dir,
				const igraph_vector_t *weights) {
  
  igraph_eit_t edgeit;
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_bool_t directed=igraph_is_directed(graph);
  igraph_vector_t degree;
  long int i;
  
  if (igraph_vector_size(weights) != no_of_edges) {
    IGRAPH_ERROR("Invalid edge weight vector length", IGRAPH_EINVAL);
  }
  
  IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, no_of_nodes));
  igraph_matrix_null(res);
  
  IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(0), &edgeit));
  IGRAPH_FINALLY(igraph_eit_destroy, &edgeit);
  
  IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
  
  if (directed) {
    IGRAPH_WARNING("Computing (weighted) Laplacian of a directed graph");

    if (dir != IGRAPH_LAPLACIAN_ROW && dir != IGRAPH_LAPLACIAN_COL) {
      IGRAPH_ERROR("Invalid direction for Laplacian", IGRAPH_EINVAL);
    }    

    if (!normalized) {

      while (!IGRAPH_EIT_END(edgeit)) {
	long int edge=IGRAPH_EIT_GET(edgeit);
	long int from=IGRAPH_FROM(graph, edge);
	long int to  =IGRAPH_TO  (graph, edge);
	igraph_real_t weight=VECTOR(*weights)[edge];
	if (weight < 0) { 
	  IGRAPH_ERROR("Negative weights are not allowed when calculating "
		       "graph Laplacian", IGRAPH_EINVAL);
	}
	if (from != to) {
	  long int wh= dir==IGRAPH_LAPLACIAN_ROW ? from : to;
	  MATRIX(*res, from, to) -= weight;
	  VECTOR(degree)[wh] += weight;
	}
	IGRAPH_EIT_NEXT(edgeit);
      }
      
      /* And the diagonal */
      for (i=0; i<no_of_nodes; i++) {
	MATRIX(*res, i, i) = VECTOR(degree)[i];
      }

    } else /* normalized */ {

      while (!IGRAPH_EIT_END(edgeit)) {
	long int edge=IGRAPH_EIT_GET(edgeit);
	long int from=IGRAPH_FROM(graph, edge);
	long int to  =IGRAPH_TO  (graph, edge);
	igraph_real_t weight=VECTOR(*weights)[edge];
	if (weight < 0) { 
	  IGRAPH_ERROR("Negative weights are not allowed when calculating "
		       "graph Laplacian", IGRAPH_EINVAL);
	}
	if (from != to) {
	  long int wh= dir==IGRAPH_LAPLACIAN_ROW ? from : to;
	  VECTOR(degree)[wh] += weight;
	}
	IGRAPH_EIT_NEXT(edgeit);
      }

      for (i=0; i<no_of_nodes; i++) {
	MATRIX(*res, i, i) = VECTOR(degree)[i] > 0 ? 1 : 0;
      }
      
      IGRAPH_EIT_RESET(edgeit);
      while (!IGRAPH_EIT_END(edgeit)) {
	long int edge=IGRAPH_EIT_GET(edgeit);
	long int from=IGRAPH_FROM(graph, edge);
	long int to  =IGRAPH_TO  (graph, edge);
	igraph_real_t weight=VECTOR(*weights)[edge];
	long int wh= dir==IGRAPH_LAPLACIAN_ROW ? from : to;
	if (from != to) {
	  MATRIX(*res, from, to) = 
	    -weight / sqrt(VECTOR(degree)[wh]);
	}
	IGRAPH_EIT_NEXT(edgeit);
      }
      
    }

  } else /* undirected */ {
    
    if (!normalized) {
      
      while (!IGRAPH_EIT_END(edgeit)) {
	long int edge=IGRAPH_EIT_GET(edgeit);
	long int from=IGRAPH_FROM(graph, edge);
	long int to  =IGRAPH_TO  (graph, edge);
	igraph_real_t weight=VECTOR(*weights)[edge];
	if (weight < 0) { 
	  IGRAPH_ERROR("Negative weights are not allowed when calculating "
		       "graph Laplacian", IGRAPH_EINVAL);
	}
	if (from != to) {
	  MATRIX(*res, from, to) -= weight;
	  MATRIX(*res, to, from) -= weight;
	  VECTOR(degree)[from] += weight;
	  VECTOR(degree)[to] += weight;
	}
	IGRAPH_EIT_NEXT(edgeit);
      }
      
      /* And the diagonal */
      for (i=0; i<no_of_nodes; i++) {
	MATRIX(*res, i, i) = VECTOR(degree)[i];
      }
      
    } else /* normalized */ {
     
      while (!IGRAPH_EIT_END(edgeit)) {
	long int edge=IGRAPH_EIT_GET(edgeit);
	long int from=IGRAPH_FROM(graph, edge);
	long int to  =IGRAPH_TO  (graph, edge);
	igraph_real_t weight=VECTOR(*weights)[edge];
	if (weight < 0) { 
	  IGRAPH_ERROR("Negative weights are not allowed when calculating "
		       "graph Laplacian", IGRAPH_EINVAL);
	}
	if (from != to) {
	  VECTOR(degree)[from] += weight;
	  VECTOR(degree)[to] += weight;
	}
	IGRAPH_EIT_NEXT(edgeit);
      }

      for (i=0; i<no_of_nodes; i++) {
	MATRIX(*res, i, i) = VECTOR(degree)[i] > 0 ? 1 : 0;
      }
      
      IGRAPH_EIT_RESET(edgeit);
      while (!IGRAPH_EIT_END(edgeit)) {
	long int edge=IGRAPH_EIT_GET(edgeit);
	long int from=IGRAPH_FROM(graph, edge);
	long int to  =IGRAPH_TO  (graph, edge);
	igraph_real_t weight=VECTOR(*weights)[edge];
	if (from != to) {
	  MATRIX(*res, from, to) = MATRIX(*res, to, from) = 
	    -weight / sqrt(VECTOR(degree)[from] * VECTOR(degree)[to]);
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

/**
 * \function igraph_laplacian
 * \brief Returns the Laplacian matrix of a graph
 * 
 * </para><para>
 * The graph Laplacian matrix is similar to an adjacency matrix but
 * contains -1's instead of 1's and the vertex degrees are included in
 * the diagonal. So the result for edge i--j is -1 if i!=j and is equal
 * to the degree of vertex i if i==j. igraph_laplacian will work on a
 * directed graph (although this does not seem to make much sense) and
 * ignores loops.
 * 
 * </para><para>
 * The normalized version of the Laplacian matrix has 1 in the diagonal and 
 * -1/sqrt(d[i]d[j]) if there is an edge from i to j.
 * 
 * </para><para>
 * If the graph is weighted (i.e. the \p weights argument is not a
 * null pointer), then all weights must be non-negative.
 * 
 * </para><para>
 * The function is supposed to work with non-simple graphs, although
 * this is not extensively tested.
 *
 * </para><para>
 * The first version of this function was written by Vincent Matossian.
 * \param graph Pointer to the graph to convert.
 * \param res Pointer to an initialized matrix object, it will be
 *        resized if needed.
 * \param normalized Whether to create a normalized Laplacian matrix.
 * \param dir Constant, it gives whether to use out degree/strength
 *        (\c IGRAPH_LAPLACIAN_ROW) or in degree/strength (\c
 *        IGRAPH_LAPLACIAN_COL) when calculating the Laplacian of a
 *        directed graph. It is ignored for undirected graphs.
 * \param weights An optional vector containing edge weights, to calculate 
 *        the weighted Laplacian matrix. Set it to a null pointer to 
 *        calculate the unweighted Laplacian.
 * \return Error code.
 *
 * Time complexity: O(|V||V|),
 * |V| is the 
 * number of vertices in the graph.
 * 
 * \sa \ref igraph_laplacian_graph(), which does the same calculation,
 * but returns the result in an igraph graph, so it might be more
 * suitable for large (and sparse) graphs.
 */

int igraph_laplacian(const igraph_t *graph, igraph_matrix_t *res,
		     igraph_bool_t normalized, 
		     igraph_laplacian_direction_t dir,
		     const igraph_vector_t *weights) {
  
  igraph_eit_t edgeit;
  long int no_of_nodes=igraph_vcount(graph);
  igraph_bool_t directed=igraph_is_directed(graph);
  long int from, to, wh;
  igraph_integer_t ffrom, fto;
  igraph_vector_t degree;  
  int i;

  if (weights) { 
    return igraph_i_weighted_laplacian(graph, res, normalized, dir, weights);
  }

  IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, no_of_nodes));
  igraph_matrix_null(res);
  IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(0), &edgeit));
  IGRAPH_FINALLY(igraph_eit_destroy, &edgeit);
  
  IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
  
  if(directed){
    igraph_neimode_t mode;

    IGRAPH_WARNING("Computing Laplacian of a directed graph");
    
    if (dir != IGRAPH_LAPLACIAN_ROW && dir != IGRAPH_LAPLACIAN_COL) {
      IGRAPH_ERROR("Invalid direction for Laplacian", IGRAPH_EINVAL);
    }
    
    mode = dir == IGRAPH_LAPLACIAN_ROW ? IGRAPH_OUT : IGRAPH_IN;

    IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(),
			       mode, IGRAPH_NO_LOOPS));
  
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
      }
      
      while (!IGRAPH_EIT_END(edgeit)) {
	igraph_edge(graph, IGRAPH_EIT_GET(edgeit), &ffrom, &fto);
	from=ffrom; to=fto;
	if (from != to) {
	  wh= dir == IGRAPH_LAPLACIAN_ROW ? from : to;
	  MATRIX(*res, from, to) = 
	    -1.0 / sqrt(VECTOR(degree)[wh]);
	}
	IGRAPH_EIT_NEXT(edgeit);
      }
    }

  } else {

    IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(),
			       IGRAPH_OUT, IGRAPH_NO_LOOPS));
  
    if (!normalized) {
      for(i=0;i<no_of_nodes;i++) {
	MATRIX(*res, i, i) = VECTOR(degree)[i];
      }
      
      while (!IGRAPH_EIT_END(edgeit)) {
	igraph_edge(graph, IGRAPH_EIT_GET(edgeit), &ffrom, &fto);
	from=ffrom;
	to=fto;	
	
	MATRIX(*res, to, from) -= 1;
	MATRIX(*res, from, to) -= 1;
	
	IGRAPH_EIT_NEXT(edgeit);
      }
    } else {
      for (i=0;i<no_of_nodes;i++) {
	MATRIX(*res, i, i) = VECTOR(degree)[i]>0 ? 1: 0;
      }
      
      while (!IGRAPH_EIT_END(edgeit)) {
	igraph_edge(graph, IGRAPH_EIT_GET(edgeit), &ffrom, &fto);
	from=ffrom; to=fto;
	if (from != to) {
	  MATRIX(*res, from, to) = MATRIX(*res, to, from) =
	    -1.0 / sqrt(VECTOR(degree)[from] * VECTOR(degree)[to]);
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

/**
 * \function igraph_laplacian_graph
 * Calculate graph Laplacian and store it in an igraph graph
 * 
 * See \ref igraph_laplacian() for more information about the graph
 * Laplacian. 
 * 
 * </para><para> 
 * Note that the current implementation does not work properly for
 * non-simple graphs, but the function does not check the input graph
 * in this respect.
 * 
 * </para><para>
 * Just like in the case of \ref igraph_laplacian(), negative weights
 * are not allowed if the graph is weighted.
 * 
 * \param graph The input graph. A warning is given if it is directed,
 *    since there is no widely accepted definition for the Laplacian
 *    of a directed graph.
 * \param res Pointer to an uninitialized graph object, the result
 *    graph will be created here.
 * \param normalized Logical, whether to calculate normalized
 *    Laplacian.
 * \param dir Constant, it gives whether to use out degree/strength
 *        (\c IGRAPH_LAPLACIAN_ROW) or in degree/strength (\c
 *        IGRAPH_LAPLACIAN_COL) when calculating the Laplacian of a
 *        directed graph. It is ignored for undirected graphs.
 * \param weights An optional weight vector for weighted graphs, or a
 *    null pointer if the graph is not weighted. Note that all weights
 *    must be non-negative.
 * \param out_weights Pointer to an initialized vector, the weights of
 *    the Laplacian graph are store here. It will be resized as
 *    needed.
 * \return Error code.
 * 
 * Time complexity: O(|V|+|E|), linear in the number of vertices plus
 * edges.
 * 
 * \sa \ref igraph_laplacian(), that performs the same calculation,
 * but returns the result in a matrix.
 */

int igraph_laplacian_graph(const igraph_t *graph, igraph_t *res,
			   igraph_bool_t normalized, 
			   igraph_laplacian_direction_t dir,
			   const igraph_vector_t *weights, 
			   igraph_vector_t *out_weights) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_vector_t el, degree;
  igraph_bool_t directed=igraph_is_directed(graph);
  long int i;
  
  IGRAPH_VECTOR_INIT_FINALLY(&el, 0);
  IGRAPH_CHECK(igraph_vector_reserve(&el, no_of_edges*2+no_of_nodes*2));
  IGRAPH_CHECK(igraph_get_edgelist(graph, &el, /*byrow=*/ 0));
  
  if (weights) {
    /* Weighted */
    
    if (igraph_vector_size(weights) != no_of_edges) {
      IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
    
    if (directed) {

      /* Weighted, directed */
      
      igraph_neimode_t mode;
      IGRAPH_WARNING("Laplacian of a directed graph");

      if (dir != IGRAPH_LAPLACIAN_ROW && dir != IGRAPH_LAPLACIAN_COL) {
	IGRAPH_ERROR("Invalid direction for Laplacian", IGRAPH_EINVAL);
      }
      
      mode = dir == IGRAPH_LAPLACIAN_ROW ? IGRAPH_OUT : IGRAPH_IN;

      IGRAPH_CHECK(igraph_strength(graph, &degree, igraph_vss_all(),
				   mode, IGRAPH_NO_LOOPS, weights));
      
      if (!normalized) {
	
	/* Weighted, directed, not normalized */
	
	IGRAPH_CHECK(igraph_vector_resize(out_weights, no_of_edges+no_of_nodes));
	for (i=0; i<no_of_edges; i++) {
	  igraph_real_t weight=VECTOR(*weights)[i];
	  if (weight < 0) { 
	    IGRAPH_ERROR("Negative weights are not allowed for calculating "
			 "the graph Laplacian", IGRAPH_EINVAL);
	  }
	  VECTOR(*out_weights)[i] = -weight;
	}
	for (i=0; i<no_of_nodes; i++) {
	  igraph_vector_push_back(&el, i);
	  igraph_vector_push_back(&el, i);
	  VECTOR(*out_weights)[no_of_edges+i] = VECTOR(degree)[i];
	}
	
      } else {

	/* Weighted, directed, normalized */

	IGRAPH_CHECK(igraph_vector_reserve(out_weights, no_of_edges+no_of_nodes));
	igraph_vector_resize(out_weights, no_of_edges);
	for (i=0; i<no_of_edges; i++) {
	  long int from=VECTOR(el)[2*i];
	  long int to=VECTOR(el)[2*i+1];
	  igraph_real_t weight=VECTOR(*weights)[i];
	  long int wh= dir==IGRAPH_LAPLACIAN_ROW ? from : to;
	  if (weight < 0) { 
	    IGRAPH_ERROR("Negative weights are not allowed for calculating "
			 "the graph Laplacian", IGRAPH_EINVAL);
	  }
	  VECTOR(*out_weights)[i] = -weight/sqrt(VECTOR(degree)[wh]);
	}
	for (i=0; i<no_of_nodes; i++) {
	  if (VECTOR(degree)[i] != 0) {
	    igraph_vector_push_back(&el, i);
	    igraph_vector_push_back(&el, i);
	    igraph_vector_push_back(out_weights, 1.0);
	  }
	}

      }

    } else {

      /* Weighted, undirected */

      IGRAPH_CHECK(igraph_strength(graph, &degree, igraph_vss_all(),
				   IGRAPH_OUT, IGRAPH_NO_LOOPS, weights));

      if (!normalized) {
	
	/* Weighted, undirected, not normalized */
	
	IGRAPH_CHECK(igraph_vector_resize(out_weights, no_of_edges+no_of_nodes));
	for (i=0; i<no_of_edges; i++) {
	  igraph_real_t weight=VECTOR(*weights)[i];
	  if (weight < 0) { 
	    IGRAPH_ERROR("Negative weights are not allowed for calculating "
			 "the graph Laplacian", IGRAPH_EINVAL);
	  }
	  VECTOR(*out_weights)[i] = -weight;
	}
	for (i=0; i<no_of_nodes; i++) {
	  igraph_vector_push_back(&el, i);
	  igraph_vector_push_back(&el, i);
	  VECTOR(*out_weights)[no_of_edges+i] = VECTOR(degree)[i];
	}
	
      } else {

	/* Weighted, undirected, normalized */
	
	IGRAPH_CHECK(igraph_vector_reserve(out_weights, no_of_edges+no_of_nodes));
	igraph_vector_resize(out_weights, no_of_edges);
	for (i=0; i<no_of_edges; i++) {
	  long int from=VECTOR(el)[2*i];
	  long int to=VECTOR(el)[2*i+1];
	  igraph_real_t weight=VECTOR(*weights)[i];
	  if (weight < 0) { 
	    IGRAPH_ERROR("Negative weights are not allowed for calculating "
			 "the graph Laplacian", IGRAPH_EINVAL);
	  }
	  VECTOR(*out_weights)[i] = 
	    -weight/sqrt(VECTOR(degree)[from] * VECTOR(degree)[to]);
	}
	for (i=0; i<no_of_nodes; i++) {
	  if (VECTOR(degree)[i] != 0) {
	    igraph_vector_push_back(&el, i);
	    igraph_vector_push_back(&el, i);
	    igraph_vector_push_back(out_weights, 1.0);
	  }
	}
      }
    }
    
  } else {

    /* Not weighted */

    IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);

    if (directed) {
      
      /* Not weighted, directed */
      
      igraph_neimode_t mode;
      IGRAPH_WARNING("Computing Laplacian of a directed graph");
      if (dir != IGRAPH_LAPLACIAN_ROW && dir != IGRAPH_LAPLACIAN_COL) {
	IGRAPH_ERROR("Invalid direction for Laplacian", IGRAPH_EINVAL);
      }
      mode = dir == IGRAPH_LAPLACIAN_ROW ? IGRAPH_OUT : IGRAPH_IN;
      IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(),
				 mode, IGRAPH_NO_LOOPS));

      if (!normalized) {

	/* Not weighted, directed, not normalized */
      
	IGRAPH_CHECK(igraph_vector_resize(out_weights, no_of_edges+no_of_nodes));
	igraph_vector_fill(out_weights, -1.0);
	for (i=0; i<no_of_nodes; i++) {
	  igraph_vector_push_back(&el, i);
	  igraph_vector_push_back(&el, i);
	  VECTOR(*out_weights)[no_of_edges+i] = VECTOR(degree)[i];
	}
      } else {

	/* Not weighted, directed, normalized */

	IGRAPH_CHECK(igraph_vector_reserve(out_weights, no_of_edges+no_of_nodes));
	igraph_vector_resize(out_weights, no_of_edges);
	for (i=0; i<no_of_edges; i++) {
	  long int idx= dir==IGRAPH_LAPLACIAN_ROW ? 2*i : 2*i+1;
	  long int wh=VECTOR(el)[idx];
	  VECTOR(*out_weights)[i] = -1.0/sqrt(VECTOR(degree)[wh]);
	}
	for (i=0; i<no_of_nodes; i++) {
	  if (VECTOR(degree)[i] != 0) {
	    igraph_vector_push_back(&el, i);
	    igraph_vector_push_back(&el, i);
	    igraph_vector_push_back(out_weights, 1.0);
	  }
	}
      }
    } else {

      /* Not weighted, undirected */

      IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(), 
				 IGRAPH_OUT, IGRAPH_NO_LOOPS));

      if (!normalized) {

	/* Not weighted, undirected, not normalized */

	IGRAPH_CHECK(igraph_vector_resize(out_weights, no_of_edges+no_of_nodes));
	igraph_vector_fill(out_weights, -1.0);
	for (i=0; i<no_of_nodes; i++) {
	  igraph_vector_push_back(&el, i);
	  igraph_vector_push_back(&el, i);
	  VECTOR(*out_weights)[no_of_edges+i] = VECTOR(degree)[i];
	}
      } else {	

	/* Not weighted, undirected, normalized */

	IGRAPH_CHECK(igraph_vector_reserve(out_weights, no_of_edges+no_of_nodes));
	igraph_vector_resize(out_weights, no_of_edges);
	for (i=0; i<no_of_edges; i++) {
	  long int from=VECTOR(el)[2*i];
	  long int to=VECTOR(el)[2*i+1];
	  VECTOR(*out_weights)[i] = 
	    -1.0/sqrt(VECTOR(degree)[from] * VECTOR(degree)[to]);
	}
	for (i=0; i<no_of_nodes; i++) {
	  if (VECTOR(degree)[i] != 0) {
	    igraph_vector_push_back(&el, i);
	    igraph_vector_push_back(&el, i);
	    igraph_vector_push_back(out_weights, 1.0);
	  }
	}
      }
    }
    
    igraph_vector_destroy(&degree);
    IGRAPH_FINALLY_CLEAN(1);
  }

  IGRAPH_CHECK(igraph_create(res, &el, no_of_nodes, directed));
  
  igraph_vector_destroy(&el);
  IGRAPH_FINALLY_CLEAN(1);
	
  return 0;
}
