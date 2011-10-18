/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2007  Gabor Csardi <csardi@rmki.kfki.hu>
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
#include "random.h"
#include "config.h"
#include <math.h>
#include <string.h>    /* memset */
#include <assert.h>

int igraph_i_eigenvector_centrality(igraph_real_t *to, const igraph_real_t *from,
				    long int n, void *extra) {
  igraph_adjlist_t *adjlist=extra;
  igraph_vector_t *neis;
  long int i, j, nlen;
  
  for (i=0; i<n; i++) {
    neis=igraph_adjlist_get(adjlist, i);
    nlen=igraph_vector_size(neis);
    to[i]=0.0;
    for (j=0; j<nlen; j++) {
      long int nei=VECTOR(*neis)[j];
      to[i] += from[nei];
    }
  }				      
  
  
  return 0;
}

typedef struct igraph_i_eigenvector_centrality_t {
  const igraph_t *graph;
  const igraph_adjedgelist_t *adjedgelist;
  const igraph_vector_t *weights;
} igraph_i_eigenvector_centrality_t;

int igraph_i_eigenvector_centrality2(igraph_real_t *to, const igraph_real_t *from,
				     long int n, void *extra) {

  igraph_i_eigenvector_centrality_t *data=extra;
  const igraph_t *graph=data->graph;
  const igraph_adjedgelist_t *adjedgelist=data->adjedgelist;
  const igraph_vector_t *weights=data->weights;
  igraph_vector_t *edges;
  long int i, j, nlen;

  for (i=0; i<n; i++) {
    edges=igraph_adjedgelist_get(adjedgelist, i);
    nlen=igraph_vector_size(edges);
    to[i]=0.0;
    for (j=0; j<nlen; j++) {
      long int edge=VECTOR(*edges)[j];
      long int nei=IGRAPH_OTHER(graph, edge, i);
      igraph_real_t w=VECTOR(*weights)[edge];
      to[i] += w * from[nei];
    }
  }

  return 0;
}

/**
 * \function igraph_eigenvector_centrality
 * Eigenvector centrality of the verices
 * 
 * Eigenvector centrality is a measure of the importance of a node in a
 * network. It assigns relative scores to all nodes in the network based
 * on the principle that connections to high-scoring nodes contribute
 * more to the score of the node in question than equal connections to
 * low-scoring nodes.
 * \param graph The input graph. It might be directed, but it will be
 *     treated as undirected anyway. 
 * \param vector Pointer to an initialized vector, it will be resized
 *     as needed. The result of the computation is stored here. It can
 *     be a null pointer, then it is ignored.
 * \param value If not a null pointer, then the eigenvalue
 *     corresponding to the found eigenvector is stored here.
 * \param scale If not zero then the result will be scaled, such that
 *     the absolute value of the maximum centrality is one.
 * \param weights A null pointer (=no edge weights), or a vector
 *     giving the weights of the edges.
 * \param options Options to ARPACK. See \ref igraph_arpack_options_t
 *    for details. Note that the function overwrites the
 *    <code>n</code> (number of vertices) parameter and 
 *    it always starts the calculation from a non-random vector
 *    calculated based on the degree of the vertices.
 * \return Error code.
 * 
 * Time complexity: depends on the input graph, usually it is O(|V|),
 * the number of vertices.
 * 
 * \sa \ref igraph_pagerank for a modification of eigenvector centrality.
 */

int igraph_eigenvector_centrality(const igraph_t *graph, igraph_vector_t *vector,
				  igraph_real_t *value, igraph_bool_t scale,
				  const igraph_vector_t *weights,
				  igraph_arpack_options_t *options) {
  
  igraph_vector_t values;
  igraph_matrix_t vectors;
  igraph_vector_t degree;
  long int i;
  
  options->n=igraph_vcount(graph);
  options->start=1;		/* no random start vector */

  if (igraph_ecount(graph) == 0) {
    /* special case: empty graph */
    if (value)
      *value = IGRAPH_NAN;
    if (vector) {
      igraph_vector_resize(vector, igraph_vcount(graph));
      igraph_vector_fill(vector, 1);
    }
	return IGRAPH_SUCCESS;
  }

  if (weights) {
	igraph_real_t min, max;

    if (igraph_vector_size(weights) != igraph_ecount(graph)) {
      IGRAPH_ERROR("Invalid length of weights vector when calculating "
                   "eigenvector centrality", IGRAPH_EINVAL);
	}
    IGRAPH_CHECK(igraph_vector_minmax(weights, &min, &max));
    if (min == 0 && max == 0) {
      /* special case: all weights are zeros */
      if (value)
        *value = IGRAPH_NAN;
      if (vector) {
        igraph_vector_resize(vector, igraph_vcount(graph));
        igraph_vector_fill(vector, 1);
      }
      return IGRAPH_SUCCESS;
    }
  }

  if (weights && igraph_is_directed(graph)) {
    IGRAPH_WARNING("Weighted directed graph in eigenvector centrality");
  }

  IGRAPH_VECTOR_INIT_FINALLY(&values, 0);
  IGRAPH_MATRIX_INIT_FINALLY(&vectors, options->n, 1);

  IGRAPH_VECTOR_INIT_FINALLY(&degree, options->n);
  IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(), 
			     IGRAPH_ALL, /*loops=*/ 0));
  for (i=0; i<options->n; i++) {
    if (VECTOR(degree)[i]) {
      MATRIX(vectors, i, 0) = VECTOR(degree)[i];
    } else {
      MATRIX(vectors, i, 0) = 1.0;
    }
  }
  igraph_vector_destroy(&degree);
  IGRAPH_FINALLY_CLEAN(1);
  
  options->n = igraph_vcount(graph);
  options->nev = 1;
  options->ncv = 3;
  options->which[0]='L'; options->which[1]='A';
  options->start=1;		/* no random start vector */

  if (!weights) {
    
    igraph_adjlist_t adjlist;

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
    
    IGRAPH_CHECK(igraph_arpack_rssolve(igraph_i_eigenvector_centrality,
				       &adjlist, options, 0, &values, &vectors));

    igraph_adjlist_destroy(&adjlist);
    IGRAPH_FINALLY_CLEAN(1);
    
  } else {
    
    igraph_adjedgelist_t adjedgelist;
    igraph_i_eigenvector_centrality_t data = { graph, &adjedgelist, weights };
    
    IGRAPH_CHECK(igraph_adjedgelist_init(graph, &adjedgelist, IGRAPH_ALL));
    IGRAPH_FINALLY(igraph_adjedgelist_destroy, &adjedgelist);
    
    IGRAPH_CHECK(igraph_arpack_rssolve(igraph_i_eigenvector_centrality2,
				       &data, options, 0, &values, &vectors));
    
    igraph_adjedgelist_destroy(&adjedgelist);
    IGRAPH_FINALLY_CLEAN(1);
  }

  if (value) {
    *value=VECTOR(values)[0];
  }
  
  if (vector) {
    igraph_real_t amax=0;
    long int which=0;
    long int i;
    IGRAPH_CHECK(igraph_vector_resize(vector, options->n));
    for (i=0; i<options->n; i++) {
      igraph_real_t tmp;
      VECTOR(*vector)[i] = MATRIX(vectors, i, 0);
      tmp=fabs(VECTOR(*vector)[i]);
      if (tmp>amax) { amax=tmp; which=i; }
    }
    if (scale && amax!=0) { igraph_vector_scale(vector, 1/VECTOR(*vector)[which]); }
  }

  if (options->info) {
    IGRAPH_WARNING("Non-zero return code from ARPACK routine!");
  }
  
  igraph_matrix_destroy(&vectors);
  igraph_vector_destroy(&values);
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}

typedef struct igraph_i_kleinberg_data_t {
  igraph_adjlist_t *in;
  igraph_adjlist_t *out;
  igraph_vector_t *tmp;
} igraph_i_kleinberg_data_t;

int igraph_i_kleinberg2(igraph_real_t *to, const igraph_real_t *from,
			long int n, void *extra) {

  igraph_adjlist_t *in = ((igraph_i_kleinberg_data_t*)extra)->in;
  igraph_adjlist_t *out = ((igraph_i_kleinberg_data_t*)extra)->out;
  igraph_vector_t *tmp = ((igraph_i_kleinberg_data_t*)extra)->tmp;
  igraph_vector_t *neis;
  long int i, j, nlen;
  
  for (i=0; i<n; i++) {
    neis=igraph_adjlist_get(in, i);
    nlen=igraph_vector_size(neis);
    VECTOR(*tmp)[i]=0.0;
    for (j=0; j<nlen; j++) {
      long int nei=VECTOR(*neis)[j];
      VECTOR(*tmp)[i] += from[nei];
    }
  }
  
  for (i=0; i<n; i++) {
    neis=igraph_adjlist_get(out, i);
    nlen=igraph_vector_size(neis);
    to[i]=0.0;
    for (j=0; j<nlen; j++) {
      long int nei=VECTOR(*neis)[j];
      to[i] += VECTOR(*tmp)[nei];
    }
  }      
  
  return 0;
}

int igraph_i_kleinberg(const igraph_t *graph, igraph_vector_t *vector,
		       igraph_real_t *value, igraph_bool_t scale,
		       igraph_arpack_options_t *options, int inout) {
  
  igraph_adjlist_t myinadjlist, myoutadjlist;
  igraph_adjlist_t *inadjlist, *outadjlist;
  igraph_vector_t tmp;
  igraph_vector_t values;
  igraph_matrix_t vectors;
  igraph_i_kleinberg_data_t extra;
  long int i;

  if (igraph_ecount(graph) == 0 || igraph_vcount(graph) == 1) {
    /* special case: empty graph or single vertex */
    if (value)
      *value = igraph_ecount(graph) ? 1.0 : IGRAPH_NAN;
    if (vector) {
      igraph_vector_resize(vector, igraph_vcount(graph));
      igraph_vector_fill(vector, 1);
    }
	return IGRAPH_SUCCESS;
  }

  options->n=igraph_vcount(graph);
  options->start=1;
  
  IGRAPH_VECTOR_INIT_FINALLY(&values, 0);
  IGRAPH_MATRIX_INIT_FINALLY(&vectors, options->n, 1);
  IGRAPH_VECTOR_INIT_FINALLY(&tmp, options->n);
  
  if (inout==0) {
    inadjlist=&myinadjlist; 
    outadjlist=&myoutadjlist;
  } else if (inout==1) {
    inadjlist=&myoutadjlist;
    outadjlist=&myinadjlist;
  } else {
    /* This should not happen */
    IGRAPH_ERROR("Invalid 'inout' argument, please do not call "
                 "this function directly", IGRAPH_FAILURE);
  }

  IGRAPH_CHECK(igraph_adjlist_init(graph, &myinadjlist, IGRAPH_IN));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &myinadjlist);
  IGRAPH_CHECK(igraph_adjlist_init(graph, &myoutadjlist, IGRAPH_OUT));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &myoutadjlist);

  IGRAPH_CHECK(igraph_degree(graph, &tmp, igraph_vss_all(), IGRAPH_ALL, 0));
  for (i=0; i<options->n; i++) {
    if (VECTOR(tmp)[i] != 0) { 
      MATRIX(vectors, i, 0) = VECTOR(tmp)[i];
    } else {
      MATRIX(vectors, i, 0) = 1.0;
    }
  }
	
  extra.in=inadjlist; extra.out=outadjlist; extra.tmp=&tmp;

  options->n = igraph_vcount(graph);
  options->nev = 1;
  options->ncv = 3;
  if (options->n < options->ncv) {
    options->ncv = options->n;
  }
  options->which[0]='L'; options->which[1]='M';
  options->start=1;		/* no random start vector */

  IGRAPH_CHECK(igraph_arpack_rssolve(igraph_i_kleinberg2, &extra,
				     options, 0, &values, &vectors));

  igraph_adjlist_destroy(&myoutadjlist);
  igraph_adjlist_destroy(&myinadjlist);
  igraph_vector_destroy(&tmp);
  IGRAPH_FINALLY_CLEAN(3);

  if (value) { 
    *value = VECTOR(values)[0];
  }

  if (vector) {
    igraph_real_t amax=0;
    long int which=0;
    long int i;
    IGRAPH_CHECK(igraph_vector_resize(vector, options->n));
    for (i=0; i<options->n; i++) {
      igraph_real_t tmp;
      VECTOR(*vector)[i] = MATRIX(vectors, i, 0);
      tmp=fabs(VECTOR(*vector)[i]);
      if (tmp>amax) { amax=tmp; which=i; }
    }
    if (scale && amax!=0) { igraph_vector_scale(vector, 1/VECTOR(*vector)[which]); }
  }
  
  if (options->info) {
    IGRAPH_WARNING("Non-zero return code from ARPACK routine!");
  }
  igraph_matrix_destroy(&vectors);
  igraph_vector_destroy(&values);
  IGRAPH_FINALLY_CLEAN(2);
  
  return 0;
}

/**
 * \function igraph_hub_score
 * Kleinberg's hub scores
 * 
 * The hub scores of the vertices are defined as the principal
 * eigenvector of <code>A*A^T</code>, where <code>A</code> is the adjacency
 * matrix of the graph, <code>A^T</code> is its transposed.
 * </para><para> 
 * See the following reference on the meaning of this score:
 * J. Kleinberg. Authoritative sources in a hyperlinked
 * environment. \emb Proc. 9th ACM-SIAM Symposium on Discrete
 * Algorithms, \eme 1998. Extended version in \emb Journal of the
 * ACM \eme 46(1999). Also appears as IBM Research Report RJ 10076, May
 * 1997.
 * \param graph The input graph. Can be directed and undirected.
 * \param vector Pointer to an initialized vector, the result is
 *    stored here. If a null pointer then it is ignored.
 * \param value If not a null pointer then the eigenvalue
 *    corresponding to the calculated eigenvector is stored here.
 * \param scale If not zero then the result will be scaled, such that
 *     the absolute value of the maximum centrality is one.
 * \param options Options to ARPACK. See \ref igraph_arpack_options_t
 *    for details. Note that the function overwrites the
 *    <code>n</code> (number of vertices) parameter and 
 *    it always starts the calculation from a non-random vector
 *    calculated based on the degree of the vertices.
 * \return Error code.
 * 
 * Time complexity: depends on the input graph, usually it is O(|V|),
 * the number of vertices.
 * 
 * \sa \ref igraph_authority_score() for the companion measure, \ref
 * igraph_pagerank(), \ref igraph_eigenvector_centrality() for similar
 * measures.
 */

int igraph_hub_score(const igraph_t *graph, igraph_vector_t *vector,
		     igraph_real_t *value, igraph_bool_t scale,
		     igraph_arpack_options_t *options) {

  return igraph_i_kleinberg(graph, vector, value, scale, options, 0);
}

/**
 * \function igraph_authority_score
 * Kleinerg's authority scores
 * 
 * The authority scores of the vertices are defined as the principal
 * eigenvector of <code>A^T*A</code>, where <code>A</code> is the adjacency
 * matrix of the graph, <code>A^T</code> is its transposed.
 * </para><para> 
 * See the following reference on the meaning of this score:
 * J. Kleinberg. Authoritative sources in a hyperlinked
 * environment. \emb Proc. 9th ACM-SIAM Symposium on Discrete
 * Algorithms, \eme 1998. Extended version in \emb Journal of the
 * ACM \eme 46(1999). Also appears as IBM Research Report RJ 10076, May
 * 1997.
 * \param graph The input graph. Can be directed and undirected.
 * \param vector Pointer to an initialized vector, the result is
 *    stored here. If a null pointer then it is ignored.
 * \param value If not a null pointer then the eigenvalue
 *    corresponding to the calculated eigenvector is stored here.
 * \param scale If not zero then the result will be scaled, such that
 *     the absolute value of the maximum centrality is one.
 * \param options Options to ARPACK. See \ref igraph_arpack_options_t
 *    for details. Note that the function overwrites the
 *    <code>n</code> (number of vertices) parameter and 
 *    it always starts the calculation from a non-random vector
 *    calculated based on the degree of the vertices.
 * \return Error code.
 * 
 * Time complexity: depends on the input graph, usually it is O(|V|),
 * the number of vertices.
 * 
 * \sa \ref igraph_hub_score() for the companion measure, \ref
 * igraph_pagerank(), \ref igraph_eigenvector_centrality() for similar
 * measures.
 */
			    
int igraph_authority_score(const igraph_t *graph, igraph_vector_t *vector,
			   igraph_real_t *value, igraph_bool_t scale,
			   igraph_arpack_options_t *options) {

  return igraph_i_kleinberg(graph, vector, value, scale, options, 1);
}

typedef struct igraph_i_pagerank_data_t {
  const igraph_t *graph;
  igraph_adjlist_t *adjlist;
  igraph_real_t damping;
  igraph_vector_t *outdegree;
  igraph_vector_t *tmp;
} igraph_i_pagerank_data_t;

typedef struct igraph_i_pagerank_data2_t {
  const igraph_t *graph;
  igraph_adjedgelist_t *adjedgelist;
  const igraph_vector_t *weights;
  igraph_real_t damping;
  igraph_vector_t *outdegree;
  igraph_vector_t *tmp;
} igraph_i_pagerank_data2_t;

int igraph_i_pagerank(igraph_real_t *to, const igraph_real_t *from,
		      long int n, void *extra) {
  
  igraph_i_pagerank_data_t *data=extra;
  igraph_adjlist_t *adjlist=data->adjlist;
  igraph_vector_t *outdegree=data->outdegree;
  igraph_vector_t *tmp=data->tmp;
  igraph_vector_t *neis;
  long int i, j, nlen;
  igraph_real_t sumfrom=0.0;

  for (i=0; i<n; i++) {
    sumfrom += from[i];
    VECTOR(*tmp)[i] = from[i] / VECTOR(*outdegree)[i];
  }
  
  for (i=0; i<n; i++) {
    neis=igraph_adjlist_get(adjlist, i);
    nlen=igraph_vector_size(neis);
    to[i]=0.0;
    for (j=0; j<nlen; j++) {
      long int nei=VECTOR(*neis)[j];
      to[i] += VECTOR(*tmp)[nei];
    }
    to[i] *= data->damping;
    to[i] += (1-data->damping)/n * sumfrom;
  }

  return 0;
}

int igraph_i_pagerank2(igraph_real_t *to, const igraph_real_t *from,
		       long int n, void *extra) {

  igraph_i_pagerank_data2_t *data=extra;
  const igraph_t *graph=data->graph;
  igraph_adjedgelist_t *adjedgelist=data->adjedgelist;
  const igraph_vector_t *weights=data->weights;
  igraph_vector_t *outdegree=data->outdegree;
  igraph_vector_t *tmp=data->tmp;
  long int i, j, nlen;
  igraph_real_t sumfrom=0.0;
  igraph_vector_t *neis;
  
  for (i=0; i<n; i++) {
    sumfrom += from[i];
    VECTOR(*tmp)[i] = from[i] / VECTOR(*outdegree)[i];
  }
  
  for (i=0; i<n; i++) {
    neis=igraph_adjedgelist_get(adjedgelist, i);
    nlen=igraph_vector_size(neis);
    to[i]=0.0;
    for (j=0; j<nlen; j++) {
      long int edge=VECTOR(*neis)[j];
      long int nei=IGRAPH_OTHER(graph, edge, i);
      to[i] += VECTOR(*weights)[edge] * VECTOR(*tmp)[nei];
    }
    to[i] *= data->damping;
    to[i] += (1-data->damping)/n * sumfrom;
  }
  
  return 0;
}

/**
 * \function igraph_pagerank
 * \brief Calculates the Google PageRank for the specified vertices.
 * 
 * This is the new PageRank implementation, based on the ARPACK
 * library. The old, power-method based implementation can be used as
 * well, it is kept under the name \ref igraph_pagerank_old().
 * 
 * </para><para>
 * Please note that the PageRank of a given vertex depends on the PageRank
 * of all other vertices, so even if you want to calculate the PageRank for
 * only some of the vertices, all of them must be calculated. Requesting
 * the PageRank for only some of the vertices does not result in any
 * performance increase at all.
 * </para>
 * <para>
 * Since the calculation is an iterative
 * process, the algorithm is stopped after a given count of iterations
 * or if the PageRank value differences between iterations are less than
 * a predefined value.
 * </para>
 * 
 * <para>
 * For the explanation of the PageRank algorithm, see the following
 * webpage:
 * http://www-db.stanford.edu/~backrub/google.html, or the
 * following reference:
 * </para>
 * 
 * <para>
 * Sergey Brin and Larry Page: The Anatomy of a Large-Scale Hypertextual
 * Web Search Engine. Proceedings of the 7th World-Wide Web Conference,
 * Brisbane, Australia, April 1998.
 * </para>
 * <para>
 * \param graph The graph object.
 * \param vector Pointer to an initialized vector, the result is
 *    stored here. It is resized as needed.
 * \param value Pointer to a real variable, the eigenvalue
 *    corresponding to the PageRank vector is stored here. It should
 *    be always exactly one.
 * \param vids The vertex ids for which the PageRank is returned.
 * \param directed Boolean, whether to consider the directedness of
 *    the edges. This is ignored for undirected graphs.
 * \param damping The damping factor ("d" in the original paper)
 * \param weights Optional edge weights, it is either a null pointer,
 *    then the edges are not weighted, or a vector of the same length
 *    as the number of edges.
 * \param options Options to ARPACK. See \ref igraph_arpack_options_t
 *    for details. Note that the function overwrites the
 *    <code>n</code> (number of vertices), <code>nev</code> (1),
 *    <code>ncv</code> (3) and <code>which</code> (LM) parameters and
 *    it always starts the calculation from a non-random vector
 *    calculated based on the degree of the vertices.
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory for
 *         temporary data. 
 *         \c IGRAPH_EINVVID, invalid vertex id in
 *         \p vids. 
 * 
 * Time complexity: TODO.
 * 
 * \sa \ref igraph_pagerank_old() for the old implementation, 
 * \ref igraph_arpack_rssolve() and \ref igraph_arpack_rnsolve() for
 * the underlying machinery.
 */

int igraph_pagerank(const igraph_t *graph, igraph_vector_t *vector,
		    igraph_real_t *value, const igraph_vs_t vids,
		    igraph_bool_t directed, igraph_real_t damping, 
		    const igraph_vector_t *weights,
		    igraph_arpack_options_t *options) {

  igraph_matrix_t values;
  igraph_matrix_t vectors;
  igraph_integer_t dirmode;
  igraph_vector_t outdegree;
  igraph_vector_t tmp;
  long int i;
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);

  if (no_of_edges == 0) {
    /* special case: empty graph */
    if (value)
      *value = IGRAPH_NAN;
    if (vector) {
      igraph_vector_resize(vector, no_of_nodes);
      igraph_vector_fill(vector, 1);
    }
	return IGRAPH_SUCCESS;
  }

  options->n = no_of_nodes;
  options->nev = 1;
  options->ncv = 3;
  options->which[0]='L'; options->which[1]='M';
  options->start=1;		/* no random start vector */

  directed = directed && igraph_is_directed(graph);

  if (weights) {
    igraph_real_t min, max;

	if (igraph_vector_size(weights) != no_of_edges) {
      IGRAPH_ERROR("Invalid length of weights vector when calculating "
                   "PageRank scores", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vector_minmax(weights, &min, &max));
    if (min == 0 && max == 0) {
      /* special case: all weights are zeros */
      if (value)
        *value = IGRAPH_NAN;
      if (vector) {
        igraph_vector_resize(vector, igraph_vcount(graph));
        igraph_vector_fill(vector, 1);
      }
      return IGRAPH_SUCCESS;
    }
  }
  
  IGRAPH_MATRIX_INIT_FINALLY(&values, 0, 0);
  IGRAPH_MATRIX_INIT_FINALLY(&vectors, options->n, 1);

  if (directed) { dirmode=IGRAPH_IN; } else { dirmode=IGRAPH_ALL; }

  IGRAPH_VECTOR_INIT_FINALLY(&outdegree, options->n);
  IGRAPH_VECTOR_INIT_FINALLY(&tmp, options->n);

  RNG_BEGIN();

  if (!weights) {
    
    igraph_adjlist_t adjlist;
    igraph_i_pagerank_data_t data = { graph, &adjlist, damping,
				      &outdegree, &tmp };

    IGRAPH_CHECK(igraph_degree(graph, &outdegree, igraph_vss_all(),
			       directed ? IGRAPH_OUT : IGRAPH_ALL, /*loops=*/ 0));
    /* Avoid division by zero */
    for (i=0; i<options->n; i++) {
      if (VECTOR(outdegree)[i]==0) {
	VECTOR(outdegree)[i]=1;
      }
      MATRIX(vectors, i, 0) = VECTOR(outdegree)[i];
    } 

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, dirmode));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
    
    IGRAPH_CHECK(igraph_arpack_rnsolve(igraph_i_pagerank,
				       &data, options, 0, &values, &vectors));

    igraph_adjlist_destroy(&adjlist);
    IGRAPH_FINALLY_CLEAN(1);
    
  } else {
    
    igraph_adjedgelist_t adjedgelist;
    igraph_i_pagerank_data2_t data = { graph, &adjedgelist, weights,
				       damping, &outdegree, &tmp };    

    IGRAPH_CHECK(igraph_adjedgelist_init(graph, &adjedgelist, dirmode));
    IGRAPH_FINALLY(igraph_adjedgelist_destroy, &adjedgelist);

    /* Weighted degree */
    for (i=0; i<no_of_edges; i++) {
      long int from=IGRAPH_FROM(graph, i);
      long int to=IGRAPH_TO(graph, i);
      igraph_real_t weight=VECTOR(*weights)[i];
      VECTOR(outdegree)[from] += weight;
      if (!directed) { 
	VECTOR(outdegree)[to]   += weight;
      }
    }
    /* Avoid division by zero */
    for (i=0; i<options->n; i++) {
      if (VECTOR(outdegree)[i]==0) {
	VECTOR(outdegree)[i]=1;
      }
      MATRIX(vectors, i, 0) = VECTOR(outdegree)[i];
    }     
    
    IGRAPH_CHECK(igraph_arpack_rnsolve(igraph_i_pagerank2,
				       &data, options, 0, &values, &vectors));
    
    igraph_adjedgelist_destroy(&adjedgelist);
    IGRAPH_FINALLY_CLEAN(1);
  }

  RNG_END();

  igraph_vector_destroy(&tmp);
  igraph_vector_destroy(&outdegree);
  IGRAPH_FINALLY_CLEAN(2);

  if (value) {
    *value=MATRIX(values, 0, 0);
  }
  
  if (vector) {
    long int i;
    igraph_vit_t vit;
    long int nodes_to_calc;
    igraph_real_t sum=0;
    
    for (i=0; i<no_of_nodes; i++) { 
      sum += MATRIX(vectors, i, 0);
    }

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    nodes_to_calc=IGRAPH_VIT_SIZE(vit);

    IGRAPH_CHECK(igraph_vector_resize(vector, nodes_to_calc));
    for (IGRAPH_VIT_RESET(vit), i=0; !IGRAPH_VIT_END(vit);
	 IGRAPH_VIT_NEXT(vit), i++) {
      VECTOR(*vector)[i] = MATRIX(vectors, (long int)IGRAPH_VIT_GET(vit), 0);
      VECTOR(*vector)[i] /= sum;
    }
    
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(1);
  }

  if (options->info) {
    IGRAPH_WARNING("Non-zero return code from ARPACK routine!");
  }
  
  igraph_matrix_destroy(&vectors);
  igraph_matrix_destroy(&values);
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}

/**
 * \ingroup structural
 * \function igraph_betweenness
 * \brief Betweenness centrality of some vertices.
 * 
 * </para><para>
 * The betweenness centrality of a vertex is the number of geodesics
 * going through it. If there are more than one geodesic between two
 * vertices, the value of these geodesics are weighted by one over the 
 * number of geodesics.
 * \param graph The graph object.
 * \param res The result of the computation, a vector containing the
 *        betweenness scores for the specified vertices.
 * \param vids The vertices of which the betweenness centrality scores
 *        will be calculated.
 * \param directed Logical, if true directed paths will be considered
 *        for directed graphs. It is ignored for undirected graphs.
 * \return Error code:
 *        \c IGRAPH_ENOMEM, not enough memory for
 *        temporary data. 
 *        \c IGRAPH_EINVVID, invalid vertex id passed in
 *        \p vids. 
 *
 * Time complexity: O(|V||E|),
 * |V| and 
 * |E| are the number of vertices and
 * edges in the graph. 
 * Note that the time complexity is independent of the number of
 * vertices for which the score is calculated.
 *
 * \sa Other centrality types: \ref igraph_degree(), \ref igraph_closeness().
 *     See \ref igraph_edge_betweenness() for calculating the betweenness score
 *     of the edges in a graph. See \ref igraph_betweenness_estimate() to
 *     estimate the betweenness score of the vertices in a graph.
 */
int igraph_betweenness(const igraph_t *graph, igraph_vector_t *res,
  const igraph_vs_t vids, igraph_bool_t directed) {
    return igraph_betweenness_estimate(graph, res, vids, directed, -1);
}

/**
 * \ingroup structural
 * \function igraph_betweenness_estimate
 * \brief Estimated betweenness centrality of some vertices.
 * 
 * </para><para>
 * The betweenness centrality of a vertex is the number of geodesics
 * going through it. If there are more than one geodesic between two
 * vertices, the value of these geodesics are weighted by one over the 
 * number of geodesics. When estimating betweenness centrality, igraph
 * takes into consideration only those paths that are shorter than or
 * equal to a prescribed length. Note that the estimated centrality
 * will always be less than the real one.
 *
 * \param graph The graph object.
 * \param res The result of the computation, a vector containing the
 *        estimated betweenness scores for the specified vertices.
 * \param vids The vertices of which the betweenness centrality scores
 *        will be estimated.
 * \param directed Logical, if true directed paths will be considered
 *        for directed graphs. It is ignored for undirected graphs.
 * \param cutoff The maximal length of paths that will be considered.
 *        If zero or negative, the exact betweenness will be calculated
 *        (no upper limit on path lengths).
 * \return Error code:
 *        \c IGRAPH_ENOMEM, not enough memory for
 *        temporary data. 
 *        \c IGRAPH_EINVVID, invalid vertex id passed in
 *        \p vids. 
 *
 * Time complexity: O(|V||E|),
 * |V| and 
 * |E| are the number of vertices and
 * edges in the graph. 
 * Note that the time complexity is independent of the number of
 * vertices for which the score is calculated.
 *
 * \sa Other centrality types: \ref igraph_degree(), \ref igraph_closeness().
 *     See \ref igraph_edge_betweenness() for calculating the betweenness score
 *     of the edges in a graph.
 */
int igraph_betweenness_estimate(const igraph_t *graph, igraph_vector_t *res, 
			const igraph_vs_t vids, igraph_bool_t directed,
                        igraph_integer_t cutoff) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  long int *distance;
  unsigned long long int *nrgeo;  /* must be long long; consider grid graphs for example */
  double *tmpscore;
  igraph_stack_t stack=IGRAPH_STACK_NULL;
  long int source;
  long int j, k;
  igraph_integer_t modein;
  igraph_vit_t vit;
  igraph_vector_t *neis;

  igraph_adjlist_t adjlist_out, adjlist_in;
  igraph_adjlist_t *adjlist_out_p, *adjlist_in_p;

  IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);

  directed=directed && igraph_is_directed(graph);
  if (directed) {
    modein=IGRAPH_IN;
    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist_out, IGRAPH_OUT));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist_out);
    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist_in, IGRAPH_IN));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist_in);
    adjlist_out_p=&adjlist_out;
    adjlist_in_p=&adjlist_in;
  } else {
    modein=IGRAPH_ALL;
    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist_out, IGRAPH_ALL));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist_out);
    adjlist_out_p=adjlist_in_p=&adjlist_out;
  }
  
  distance=igraph_Calloc(no_of_nodes, long int);
  if (distance==0) {
    IGRAPH_ERROR("betweenness failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, distance);
  nrgeo=igraph_Calloc(no_of_nodes, unsigned long long int);
  if (nrgeo==0) {
    IGRAPH_ERROR("betweenness failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, nrgeo);
  tmpscore=igraph_Calloc(no_of_nodes, double);
  if (tmpscore==0) {
    IGRAPH_ERROR("betweenness failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, tmpscore);

  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
  igraph_stack_init(&stack, no_of_nodes);
  IGRAPH_FINALLY(igraph_stack_destroy, &stack);
    
  IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_VIT_SIZE(vit)));
  igraph_vector_null(res);

  /* here we go */
  
  for (source=0; source<no_of_nodes; source++) {
    IGRAPH_PROGRESS("Betweenness centrality: ", 100.0*source/no_of_nodes, 0);
    IGRAPH_ALLOW_INTERRUPTION();

    memset(distance, 0, no_of_nodes*sizeof(long int));
    memset(nrgeo, 0, no_of_nodes*sizeof(unsigned long long int));
    memset(tmpscore, 0, no_of_nodes*sizeof(double));
    igraph_stack_clear(&stack); /* it should be empty anyway... */
    
    IGRAPH_CHECK(igraph_dqueue_push(&q, source));
    nrgeo[source]=1;
    distance[source]=0;
    
    while (!igraph_dqueue_empty(&q)) {
      long int actnode=igraph_dqueue_pop(&q);

      if (cutoff > 0 && distance[actnode] >= cutoff) continue;
       
      neis = igraph_adjlist_get(adjlist_out_p, actnode);
      for (j=0; j<igraph_vector_size(neis); j++) {
        long int neighbor=VECTOR(*neis)[j];
        if (nrgeo[neighbor] != 0) {
	      /* we've already seen this node, another shortest path? */
	      if (distance[neighbor]==distance[actnode]+1) {
	        nrgeo[neighbor]+=nrgeo[actnode];
	      }
	    } else {
	      /* we haven't seen this node yet */
	      nrgeo[neighbor]+=nrgeo[actnode];
              distance[neighbor]=distance[actnode]+1;
	      IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
	      IGRAPH_CHECK(igraph_stack_push(&stack, neighbor));
	    }
      }
    } /* while !igraph_dqueue_empty */

    /* Ok, we've the distance of each node and also the number of
       shortest paths to them. Now we do an inverse search, starting
       with the farthest nodes. */
    while (!igraph_stack_empty(&stack)) {
      long int actnode=igraph_stack_pop(&stack);      
      if (distance[actnode]<=1) { continue; } /* skip source node */
      
      /* set the temporary score of the friends */
      neis = igraph_adjlist_get(adjlist_in_p, actnode);
      for (j=0; j<igraph_vector_size(neis); j++) {
        long int neighbor=VECTOR(*neis)[j];
	    if (distance[neighbor]==distance[actnode]-1 && nrgeo[neighbor] != 0) {
	      tmpscore[neighbor] += 
	        (tmpscore[actnode]+1)*nrgeo[neighbor]/nrgeo[actnode];
	    }
      }
    }
    
    /* Ok, we've the scores for this source */
    for (k=0, IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); 
	 IGRAPH_VIT_NEXT(vit), k++) {
      long int node=IGRAPH_VIT_GET(vit);
      VECTOR(*res)[k] += tmpscore[node];
      tmpscore[node] = 0.0; /* in case a node is in vids multiple times */
    }

  } /* for source < no_of_nodes */

  /* divide by 2 for undirected graph */
  if (!directed) {
    for (j=0; j<igraph_vector_size(res); j++) {
      VECTOR(*res)[j] /= 2.0;
    }
  }
  
  /* clean  */
  igraph_Free(distance);
  igraph_Free(nrgeo);
  igraph_Free(tmpscore);
  
  igraph_dqueue_destroy(&q);
  igraph_stack_destroy(&stack);
  igraph_vit_destroy(&vit);
  IGRAPH_FINALLY_CLEAN(6);

  if (directed) {
    igraph_adjlist_destroy(&adjlist_out);
    igraph_adjlist_destroy(&adjlist_in);
    IGRAPH_FINALLY_CLEAN(2);
  } else {
    igraph_adjlist_destroy(&adjlist_out);
    IGRAPH_FINALLY_CLEAN(1);
  }

  return 0;
}

/**
 * \ingroup structural
 * \function igraph_edge_betweenness
 * \brief Betweenness centrality of the edges.
 * 
 * </para><para>
 * The betweenness centrality of an edge is the number of geodesics
 * going through it. If there are more than one geodesics between two
 * vertices, the value of these geodesics are weighted by one over the 
 * number of geodesics.
 * \param graph The graph object.
 * \param result The result of the computation, vector containing the
 *        betweenness scores for the edges.
 * \param directed Logical, if true directed paths will be considered
 *        for directed graphs. It is ignored for undirected graphs.
 * \return Error code:
 *        \c IGRAPH_ENOMEM, not enough memory for
 *        temporary data. 
 *
 * Time complexity: O(|V||E|),
 * |V| and
 * |E| are the number of vertices and
 * edges in the graph. 
 *
 * \sa Other centrality types: \ref igraph_degree(), \ref igraph_closeness().
 *     See \ref igraph_edge_betweenness() for calculating the betweenness score
 *     of the edges in a graph. See \ref igraph_edge_betweenness_estimate() to
 *     estimate the betweenness score of the edges in a graph.
 */
int igraph_edge_betweenness(const igraph_t *graph, igraph_vector_t *result,
                            igraph_bool_t directed) {
  return igraph_edge_betweenness_estimate(graph, result, directed, -1);
}

/**
 * \ingroup structural
 * \function igraph_edge_betweenness_estimate
 * \brief Estimated betweenness centrality of the edges.
 * 
 * </para><para>
 * The betweenness centrality of an edge is the number of geodesics
 * going through it. If there are more than one geodesics between two
 * vertices, the value of these geodesics are weighted by one over the 
 * number of geodesics. When estimating betweenness centrality, igraph
 * takes into consideration only those paths that are shorter than or
 * equal to a prescribed length. Note that the estimated centrality
 * will always be less than the real one.
 * \param graph The graph object.
 * \param result The result of the computation, vector containing the
 *        betweenness scores for the edges.
 * \param directed Logical, if true directed paths will be considered
 *        for directed graphs. It is ignored for undirected graphs.
 * \param cutoff The maximal length of paths that will be considered.
 *        If zero or negative, the exact betweenness will be calculated
 *        (no upper limit on path lengths).
 * \return Error code:
 *        \c IGRAPH_ENOMEM, not enough memory for
 *        temporary data. 
 *
 * Time complexity: O(|V||E|),
 * |V| and
 * |E| are the number of vertices and
 * edges in the graph. 
 *
 * \sa Other centrality types: \ref igraph_degree(), \ref igraph_closeness().
 *     See \ref igraph_betweenness() for calculating the betweenness score
 *     of the vertices in a graph.
 */
int igraph_edge_betweenness_estimate(const igraph_t *graph, igraph_vector_t *result,
                                     igraph_bool_t directed, igraph_integer_t cutoff) {
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  long int *distance;
  unsigned long long int *nrgeo;
  double *tmpscore;
  igraph_stack_t stack=IGRAPH_STACK_NULL;
  long int source;
  long int j;

  igraph_adjedgelist_t elist_out, elist_in;
  igraph_adjedgelist_t *elist_out_p, *elist_in_p;
  igraph_vector_t *neip;
  long int neino;
  long int i;
  igraph_integer_t modein;

  directed=directed && igraph_is_directed(graph);
  if (directed) {
    modein=IGRAPH_IN;
    IGRAPH_CHECK(igraph_adjedgelist_init(graph, &elist_out, IGRAPH_OUT));
    IGRAPH_FINALLY(igraph_adjedgelist_destroy, &elist_out);
    IGRAPH_CHECK(igraph_adjedgelist_init(graph, &elist_in, IGRAPH_IN));
    IGRAPH_FINALLY(igraph_adjedgelist_destroy, &elist_in);
    elist_out_p=&elist_out;
    elist_in_p=&elist_in;
  } else {
    modein=IGRAPH_ALL;
    IGRAPH_CHECK(igraph_adjedgelist_init(graph,&elist_out, IGRAPH_ALL));
    IGRAPH_FINALLY(igraph_adjedgelist_destroy, &elist_out);
    elist_out_p=elist_in_p=&elist_out;
  }
  
  distance=igraph_Calloc(no_of_nodes, long int);
  if (distance==0) {
    IGRAPH_ERROR("edge betweenness failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, distance);
  nrgeo=igraph_Calloc(no_of_nodes, unsigned long long int);
  if (nrgeo==0) {
    IGRAPH_ERROR("edge betweenness failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, nrgeo);
  tmpscore=igraph_Calloc(no_of_nodes, double);
  if (tmpscore==0) {
    IGRAPH_ERROR("edge betweenness failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, tmpscore);

  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
  IGRAPH_CHECK(igraph_stack_init(&stack, no_of_nodes));
  IGRAPH_FINALLY(igraph_stack_destroy, &stack);

  IGRAPH_CHECK(igraph_vector_resize(result, no_of_edges));

  igraph_vector_null(result);

  /* here we go */
  
  for (source=0; source<no_of_nodes; source++) {
    IGRAPH_PROGRESS("Edge betweenness centrality: ", 100.0*source/no_of_nodes, 0);
    IGRAPH_ALLOW_INTERRUPTION();

    memset(distance, 0, no_of_nodes*sizeof(long int));
    memset(nrgeo, 0, no_of_nodes*sizeof(unsigned long long int));
    memset(tmpscore, 0, no_of_nodes*sizeof(double));
    igraph_stack_clear(&stack); /* it should be empty anyway... */
    
    IGRAPH_CHECK(igraph_dqueue_push(&q, source));
      
    nrgeo[source]=1;
    distance[source]=0;
    
    while (!igraph_dqueue_empty(&q)) {
      long int actnode=igraph_dqueue_pop(&q);

      if (cutoff > 0 && distance[actnode] >= cutoff ) continue;

      neip=igraph_adjedgelist_get(elist_out_p, actnode);
      neino=igraph_vector_size(neip);
      for (i=0; i<neino; i++) {
	igraph_integer_t edge=VECTOR(*neip)[i], from, to;
	long int neighbor;
	igraph_edge(graph, edge, &from, &to);
	neighbor = actnode!=from ? from : to;
	if (nrgeo[neighbor] != 0) {
	  /* we've already seen this node, another shortest path? */
	  if (distance[neighbor]==distance[actnode]+1) {
	    nrgeo[neighbor]+=nrgeo[actnode];
	  }
	} else {
	  /* we haven't seen this node yet */
	  nrgeo[neighbor]+=nrgeo[actnode];
	  distance[neighbor]=distance[actnode]+1;
	  IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
	  IGRAPH_CHECK(igraph_stack_push(&stack, neighbor));
	}
      }
    } /* while !igraph_dqueue_empty */
    
    /* Ok, we've the distance of each node and also the number of
       shortest paths to them. Now we do an inverse search, starting
       with the farthest nodes. */
    while (!igraph_stack_empty(&stack)) {
      long int actnode=igraph_stack_pop(&stack);
      if (distance[actnode]<1) { continue; } /* skip source node */
      
      /* set the temporary score of the friends */
      neip=igraph_adjedgelist_get(elist_in_p, actnode);
      neino=igraph_vector_size(neip);
      for (i=0; i<neino; i++) {
	igraph_integer_t from, to;
	long int neighbor;
	long int edgeno=VECTOR(*neip)[i];
	igraph_edge(graph, edgeno, &from, &to);
	neighbor= actnode != from ? from : to;
	if (distance[neighbor]==distance[actnode]-1 &&
	    nrgeo[neighbor] != 0) {
	  tmpscore[neighbor] +=
	    (tmpscore[actnode]+1)*nrgeo[neighbor]/nrgeo[actnode];
	  VECTOR(*result)[edgeno] +=
	    (tmpscore[actnode]+1)*nrgeo[neighbor]/nrgeo[actnode];
	}
      }
    }
    /* Ok, we've the scores for this source */
  } /* for source <= no_of_nodes */
  IGRAPH_PROGRESS("Edge betweenness centrality: ", 100.0, 0);

  /* clean and return */
  igraph_Free(distance);
  igraph_Free(nrgeo);
  igraph_Free(tmpscore);
  igraph_dqueue_destroy(&q);
  igraph_stack_destroy(&stack);
  IGRAPH_FINALLY_CLEAN(5);

  if (directed) {
    igraph_adjedgelist_destroy(&elist_out);
    igraph_adjedgelist_destroy(&elist_in);
    IGRAPH_FINALLY_CLEAN(2);
  } else {
    igraph_adjedgelist_destroy(&elist_out);
    IGRAPH_FINALLY_CLEAN(1);
  }

  /* divide by 2 for undirected graph */
  if (!directed || !igraph_is_directed(graph)) {
    for (j=0; j<igraph_vector_size(result); j++) {
      VECTOR(*result)[j] /= 2.0;
    }
  }
  
  return 0;
}

/**
 * \ingroup structural
 * \function igraph_closeness
 * \brief Closeness centrality calculations for some vertices.
 *
 * </para><para>
 * The closeness centrality of a vertex measures how easily other
 * vertices can be reached from it (or the other way: how easily it
 * can be reached from the other vertices). It is defined as the
 * number of the number of vertices minus one divided by the sum of the
 * lengths of all geodesics from/to the given vertex.
 *
 * </para><para>
 * If the graph is not connected, and there is no path between two
 * vertices, the number of vertices is used instead the length of the
 * geodesic. This is always longer than the longest possible geodesic.
 * 
 * \param graph The graph object.
 * \param res The result of the computation, a vector containing the
 *        closeness centrality scores for the given vertices.
 * \param vids Vector giving the vertices for which the closeness
 *        centrality scores will be computed.
 * \param mode The type of shortest paths to be used for the
 *        calculation in directed graphs. Possible values: 
 *        \clist
 *        \cli IGRAPH_OUT 
 *          the lengths of the outgoing paths are calculated. 
 *        \cli IGRAPH_IN 
 *          the lengths of the incoming paths are calculated. 
 *        \cli IGRAPH_ALL
 *          the directed graph is considered as an
 *          undirected one for the computation.
 *        \endclist
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *           not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           invalid vertex id passed.
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 *
 * Time complexity: O(n|E|),
 * n is the number 
 * of vertices for which the calculation is done and
 * |E| is the number 
 * of edges in the graph.
 *
 * \sa Other centrality types: \ref igraph_degree(), \ref igraph_betweenness().
 *   See \ref igraph_closeness_estimate() to estimate closeness values.
 */
int igraph_closeness(const igraph_t *graph, igraph_vector_t *res,
                     const igraph_vs_t vids, igraph_neimode_t mode) {
  return igraph_closeness_estimate(graph, res, vids, mode, -1);
}

/**
 * \ingroup structural
 * \function igraph_closeness_estimate
 * \brief Closeness centrality estimations for some vertices.
 *
 * </para><para>
 * The closeness centrality of a vertex measures how easily other
 * vertices can be reached from it (or the other way: how easily it
 * can be reached from the other vertices). It is defined as the
 * number of the number of vertices minus one divided by the sum of the
 * lengths of all geodesics from/to the given vertex. When estimating
 * closeness centrality, igraph considers paths having a length less than
 * or equal to a prescribed cutoff value.
 *
 * </para><para>
 * If the graph is not connected, and there is no such path between two
 * vertices, the number of vertices is used instead the length of the
 * geodesic. This is always longer than the longest possible geodesic.
 *
 * </para><para>
 * Since the estimation considers vertex pairs with a distance greater than
 * the given value as disconnected, the resulting estimation will always be
 * lower than the actual closeness centrality.
 * 
 * \param graph The graph object.
 * \param res The result of the computation, a vector containing the
 *        closeness centrality scores for the given vertices.
 * \param vids Vector giving the vertices for which the closeness
 *        centrality scores will be computed.
 * \param mode The type of shortest paths to be used for the
 *        calculation in directed graphs. Possible values: 
 *        \clist
 *        \cli IGRAPH_OUT 
 *          the lengths of the outgoing paths are calculated. 
 *        \cli IGRAPH_IN 
 *          the lengths of the incoming paths are calculated. 
 *        \cli IGRAPH_ALL
 *          the directed graph is considered as an
 *          undirected one for the computation.
 *        \endclist
 * \param cutoff The maximal length of paths that will be considered.
 *        If zero or negative, the exact closeness will be calculated
 *        (no upper limit on path lengths).
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *           not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           invalid vertex id passed.
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 *
 * Time complexity: O(n|E|),
 * n is the number 
 * of vertices for which the calculation is done and
 * |E| is the number 
 * of edges in the graph.
 *
 * \sa Other centrality types: \ref igraph_degree(), \ref igraph_betweenness().
 */
int igraph_closeness_estimate(const igraph_t *graph, igraph_vector_t *res, 
		              const igraph_vs_t vids, igraph_neimode_t mode,
                              igraph_integer_t cutoff) {
  long int no_of_nodes=igraph_vcount(graph);
  igraph_vector_t already_counted, *neis;
  long int i, j;
  long int nodes_reached;
  igraph_adjlist_t allneis;

  igraph_dqueue_t q;
  
  long int nodes_to_calc;
  igraph_vit_t vit;

  IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);

  nodes_to_calc=IGRAPH_VIT_SIZE(vit);
  
  if (mode != IGRAPH_OUT && mode != IGRAPH_IN && 
      mode != IGRAPH_ALL) {
    IGRAPH_ERROR("calculating closeness", IGRAPH_EINVMODE);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&already_counted, no_of_nodes);
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);

  IGRAPH_CHECK(igraph_adjlist_init(graph, &allneis, mode));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &allneis);

  IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));
  igraph_vector_null(res);
  
  for (IGRAPH_VIT_RESET(vit), i=0; 
       !IGRAPH_VIT_END(vit); 
       IGRAPH_VIT_NEXT(vit), i++) {
    IGRAPH_CHECK(igraph_dqueue_push(&q, IGRAPH_VIT_GET(vit)));
    IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
    nodes_reached=1;
    VECTOR(already_counted)[(long int)IGRAPH_VIT_GET(vit)]=i+1;

    IGRAPH_PROGRESS("Closeness: ", 100.0*i/no_of_nodes, NULL);
    IGRAPH_ALLOW_INTERRUPTION();
    
    while (!igraph_dqueue_empty(&q)) {
      long int act=igraph_dqueue_pop(&q);
      long int actdist=igraph_dqueue_pop(&q);
      
      VECTOR(*res)[i] += actdist;

      if (cutoff>0 && actdist>=cutoff) continue;

      neis=igraph_adjlist_get(&allneis, act);
      for (j=0; j<igraph_vector_size(neis); j++) {
        long int neighbor=VECTOR(*neis)[j];
        if (VECTOR(already_counted)[neighbor] == i+1) { continue; }
        VECTOR(already_counted)[neighbor] = i+1;
        nodes_reached++;
        IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
        IGRAPH_CHECK(igraph_dqueue_push(&q, actdist+1));
      }
    }
    VECTOR(*res)[i] += ((igraph_integer_t)no_of_nodes * (no_of_nodes-nodes_reached));
    VECTOR(*res)[i] = (no_of_nodes-1) / VECTOR(*res)[i];
  }

  IGRAPH_PROGRESS("Closeness: ", 100.0, NULL);

  /* Clean */
  igraph_dqueue_destroy(&q);
  igraph_vector_destroy(&already_counted);
  igraph_vit_destroy(&vit);
  igraph_adjlist_destroy(&allneis);
  IGRAPH_FINALLY_CLEAN(4);
  
  return 0;
}


