/* -*- mode: C -*-  */
/* vim:set ts=2 sts=2 sw=2 et: */
/* 
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include <math.h>
#include <string.h>    /* memset */
#include <assert.h>
#include "igraph_centrality.h"
#include "igraph_math.h"
#include "igraph_memory.h"
#include "igraph_random.h"
#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_progress.h"
#include "igraph_interrupt_internal.h"
#include "igraph_topology.h"
#include "igraph_types_internal.h"
#include "igraph_stack.h"
#include "igraph_dqueue.h"
#include "config.h"

#include "bigint.h"
#include "prpack.h"

int igraph_personalized_pagerank_arpack(const igraph_t *graph, 
		    igraph_vector_t *vector,
		    igraph_real_t *value, const igraph_vs_t vids,
		    igraph_bool_t directed, igraph_real_t damping, 
		    igraph_vector_t *reset,
		    const igraph_vector_t *weights,
		    igraph_arpack_options_t *options);

igraph_bool_t igraph_i_vector_mostly_negative(const igraph_vector_t *vector) {
  /* Many of the centrality measures correspond to the eigenvector of some
   * matrix. When v is an eigenvector, c*v is also an eigenvector, therefore
   * it may happen that all the scores in the eigenvector are negative, in which
   * case we want to negate them since the centrality scores should be positive.
   * However, since ARPACK is not always stable, sometimes it happens that
   * *some* of the centrality scores are small negative numbers. This function
   * helps distinguish between the two cases; it should return true if most of
   * the values are relatively large negative numbers, in which case we should
   * negate the eigenvector.
   */
  long int i, n = igraph_vector_size(vector);
  igraph_real_t mi, ma;

  if (n == 0)
    return 0;

  mi = ma = VECTOR(*vector)[0];
  for (i = 1; i < n; i++) {
    if (VECTOR(*vector)[i] < mi)
      mi = VECTOR(*vector)[i];
    if (VECTOR(*vector)[i] > ma)
      ma = VECTOR(*vector)[i];
  }

  if (mi >= 0)
    return 0;
  if (ma <= 0)
    return 1;

  mi /= ma;
  return (mi < 1e-5) ? 1 : 0;
}

int igraph_i_eigenvector_centrality(igraph_real_t *to, const igraph_real_t *from,
				    int n, void *extra) {
  igraph_adjlist_t *adjlist=extra;
  igraph_vector_int_t *neis;
  long int i, j, nlen;
  
  for (i=0; i<n; i++) {
    neis=igraph_adjlist_get(adjlist, i);
    nlen=igraph_vector_int_size(neis);
    to[i]=0.0;
    for (j=0; j<nlen; j++) {
      long int nei=(long int) VECTOR(*neis)[j];
      to[i] += from[nei];
    }
  }				      
  
  
  return 0;
}

typedef struct igraph_i_eigenvector_centrality_t {
  const igraph_t *graph;
  const igraph_inclist_t *inclist;
  const igraph_vector_t *weights;
} igraph_i_eigenvector_centrality_t;

int igraph_i_eigenvector_centrality2(igraph_real_t *to, const igraph_real_t *from,
				     int n, void *extra) {

  igraph_i_eigenvector_centrality_t *data=extra;
  const igraph_t *graph=data->graph;
  const igraph_inclist_t *inclist=data->inclist;
  const igraph_vector_t *weights=data->weights;
  igraph_vector_int_t *edges;
  long int i, j, nlen;

  for (i=0; i<n; i++) {
    edges=igraph_inclist_get(inclist, i);
    nlen=igraph_vector_int_size(edges);
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

int igraph_i_eigenvector_centrality_loop(igraph_adjlist_t *adjlist) {
  
  long int i, j, k, nlen, n=igraph_adjlist_size(adjlist);
  igraph_vector_int_t *neis;

  for (i=0; i<n; i++) {
    neis=igraph_adjlist_get(adjlist, i);
    nlen=igraph_vector_int_size(neis);    
    for (j=0; j<nlen && VECTOR(*neis)[j]<i; j++) ;
    for (k=j; k<nlen && VECTOR(*neis)[k]==i; k++) ;
    if (k!=j) {
      /* First loop edge is 'j', first non-loop edge is 'k' */
      igraph_vector_int_remove_section(neis, j+(k-j)/2, k);
    }
  }

  return 0;
}

int igraph_eigenvector_centrality_undirected(const igraph_t *graph, igraph_vector_t *vector,
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
      *value = 0;
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
        *value = 0;
      if (vector) {
        igraph_vector_resize(vector, igraph_vcount(graph));
        igraph_vector_fill(vector, 1);
      }
      return IGRAPH_SUCCESS;
    }
  }

  IGRAPH_VECTOR_INIT_FINALLY(&values, 0);
  IGRAPH_MATRIX_INIT_FINALLY(&vectors, options->n, 1);

  IGRAPH_VECTOR_INIT_FINALLY(&degree, options->n);
  IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(), 
			     IGRAPH_ALL, /*loops=*/ 0));
  RNG_BEGIN();
  for (i=0; i<options->n; i++) {
    if (VECTOR(degree)[i]) {
      MATRIX(vectors, i, 0) = VECTOR(degree)[i] + RNG_UNIF(-1e-4, 1e-4);
    } else {
      MATRIX(vectors, i, 0) = 1.0;
    }
  }
  RNG_END();
  igraph_vector_destroy(&degree);
  IGRAPH_FINALLY_CLEAN(1);
  
  options->n = igraph_vcount(graph);
  options->nev = 1;
  options->ncv = 0;   /* 0 means "automatic" in igraph_arpack_rssolve */
  options->which[0]='L'; options->which[1]='A';
  options->start=1;		/* no random start vector */

  if (!weights) {
    
    igraph_adjlist_t adjlist;

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    IGRAPH_CHECK(igraph_i_eigenvector_centrality_loop(&adjlist));
    
    IGRAPH_CHECK(igraph_arpack_rssolve(igraph_i_eigenvector_centrality,
				       &adjlist, options, 0, &values, &vectors));

    igraph_adjlist_destroy(&adjlist);
    IGRAPH_FINALLY_CLEAN(1);
    
  } else {
    
    igraph_inclist_t inclist;
    igraph_i_eigenvector_centrality_t data = { graph, &inclist, weights };
    
    IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, IGRAPH_ALL));
    IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);

    IGRAPH_CHECK(igraph_inclist_remove_duplicate(graph, &inclist));
    
    IGRAPH_CHECK(igraph_arpack_rssolve(igraph_i_eigenvector_centrality2,
				       &data, options, 0, &values, &vectors));
    
    igraph_inclist_destroy(&inclist);
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

    if (VECTOR(values)[0] <= 0) {
      /* Pathological case: largest eigenvalue is zero, therefore all the
       * scores can also be zeros, this will be a valid eigenvector.
       * This usually happens with graphs that have lots of sinks and
       * sources only. */
      igraph_vector_fill(vector, 0);
    } else {
      for (i=0; i<options->n; i++) {
        igraph_real_t tmp;
        VECTOR(*vector)[i] = MATRIX(vectors, i, 0);
        tmp=fabs(VECTOR(*vector)[i]);
        if (tmp>amax) { amax=tmp; which=i; }
      }
      if (scale && amax!=0) { 
        igraph_vector_scale(vector, 1/VECTOR(*vector)[which]); 
      } else if (igraph_i_vector_mostly_negative(vector)) {
        igraph_vector_scale(vector, -1.0);
      }

      /* Correction for numeric inaccuracies (eliminating -0.0) */
      for (i=0; i<options->n; i++) {
        if (VECTOR(*vector)[i] < 0)
		    VECTOR(*vector)[i] = 0;
      }
    }
  }

  if (options->info) {
    IGRAPH_WARNING("Non-zero return code from ARPACK routine!");
  }
  
  igraph_matrix_destroy(&vectors);
  igraph_vector_destroy(&values);
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}

/* int igraph_i_evcent_dir(igraph_real_t *to, const igraph_real_t *from, */
/* 			long int n, void *extra) { */
/*   /\* TODO *\/ */
/*   return 0; */
/* } */

/* int igraph_i_evcent_dir2(igraph_real_t *to, const igraph_real_t *from, */
/* 			 long int n, void *extra) { */
/*   /\* TODO *\/ */
/*   return 0; */
/* } */

int igraph_eigenvector_centrality_directed(const igraph_t *graph, igraph_vector_t *vector,
					   igraph_real_t *value, igraph_bool_t scale,
					   const igraph_vector_t *weights,
					   igraph_arpack_options_t *options) {
  
  igraph_matrix_t values;
  igraph_matrix_t vectors;
  igraph_vector_t indegree;
  igraph_bool_t dag;
  long int i;

  if (igraph_ecount(graph) == 0) {
    /* special case: empty graph */
    if (value)
      *value = 0;
    if (vector) {
      igraph_vector_resize(vector, igraph_vcount(graph));
      igraph_vector_fill(vector, 1);
    }
    return IGRAPH_SUCCESS;
  }

  /* Quick check: if the graph is a DAG, all the eigenvector centralities are
   * zeros, and so is the eigenvalue */
  IGRAPH_CHECK(igraph_is_dag(graph, &dag));
  if (dag) {
    /* special case: graph is a DAG */
    IGRAPH_WARNING("graph is directed and acyclic; eigenvector centralities "
        "will be zeros");
    if (value)
      *value = 0;
    if (vector) {
      igraph_vector_resize(vector, igraph_vcount(graph));
      igraph_vector_fill(vector, 0);
    }
    return IGRAPH_SUCCESS;
  }

  if (weights) {
    igraph_real_t min, max;

    if (igraph_vector_size(weights) != igraph_ecount(graph)) {
      IGRAPH_ERROR("Invalid length of weights vector when calculating "
                   "eigenvector centrality", IGRAPH_EINVAL);
    }
    if (igraph_is_directed(graph)) {
      IGRAPH_WARNING("Weighted directed graph in eigenvector centrality");
    }

	IGRAPH_CHECK(igraph_vector_minmax(weights, &min, &max));

	if (min < 0.0) {
      IGRAPH_WARNING("Negative weights, eigenpair might be complex");
    }
    if (min == 0.0 && max == 0.0) {
      /* special case: all weights are zeros */
      if (value)
        *value = 0;
      if (vector) {
        igraph_vector_resize(vector, igraph_vcount(graph));
        igraph_vector_fill(vector, 1);
      }
      return IGRAPH_SUCCESS;
    }
  }

  options->n=igraph_vcount(graph);
  options->start=1;
  options->nev=1;
  options->ncv = 0;   /* 0 means "automatic" in igraph_arpack_rnsolve */
  /* LM mode is not OK here because +1 and -1 can be eigenvalues at the
   * same time, e.g.: a -> b -> a, c -> a */
  options->which[0]='L' ; options->which[1]='R';

  IGRAPH_MATRIX_INIT_FINALLY(&values, 0, 0);
  IGRAPH_MATRIX_INIT_FINALLY(&vectors, options->n, 1);
  
  IGRAPH_VECTOR_INIT_FINALLY(&indegree, options->n);
  IGRAPH_CHECK(igraph_strength(graph, &indegree, igraph_vss_all(), 
	       IGRAPH_IN, /*loops=*/ 1, weights));
  RNG_BEGIN();
  for (i=0; i<options->n; i++) {
    if (VECTOR(indegree)[i]) {
      MATRIX(vectors, i, 0) = VECTOR(indegree)[i] + RNG_UNIF(-1e-4, 1e-4);
    } else {
      MATRIX(vectors, i, 0) = 1.0;
    }
  }
  RNG_END();
  igraph_vector_destroy(&indegree);
  IGRAPH_FINALLY_CLEAN(1);

  if (!weights) {
    igraph_adjlist_t adjlist;

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_IN));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
    
    IGRAPH_CHECK(igraph_arpack_rnsolve(igraph_i_eigenvector_centrality,
				       &adjlist, options, 0, &values, 
				       &vectors));
    
    igraph_adjlist_destroy(&adjlist);
    IGRAPH_FINALLY_CLEAN(1);
  } else {
    igraph_inclist_t inclist;
    igraph_i_eigenvector_centrality_t data={ graph, &inclist, weights };

    IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, IGRAPH_IN));
    IGRAPH_FINALLY(igraph_inclist_destroy, &inclist); 
    
    IGRAPH_CHECK(igraph_arpack_rnsolve(igraph_i_eigenvector_centrality2,
				       &data, options, 0, &values, &vectors));
    
    igraph_inclist_destroy(&inclist);
    IGRAPH_FINALLY_CLEAN(1);
  }

  if (value) {
    *value=MATRIX(values, 0, 0);
  }

  if (vector) {
    igraph_real_t amax=0;
    long int which=0;
    long int i;
    IGRAPH_CHECK(igraph_vector_resize(vector, options->n));

    if (MATRIX(values, 0, 0) <= 0) {
      /* Pathological case: largest eigenvalue is zero, therefore all the
       * scores can also be zeros, this will be a valid eigenvector.
       * This usually happens with graphs that have lots of sinks and
       * sources only. */
      igraph_vector_fill(vector, 0);
      MATRIX(values, 0, 0) = 0;
    } else {
      for (i=0; i<options->n; i++) {
        igraph_real_t tmp;
        VECTOR(*vector)[i] = MATRIX(vectors, i, 0);
        tmp=fabs(VECTOR(*vector)[i]);
        if (tmp>amax) { amax=tmp; which=i; }
      }
      if (scale && amax!=0) { 
        igraph_vector_scale(vector, 1/VECTOR(*vector)[which]); 
      } else if (igraph_i_vector_mostly_negative(vector)) {
        igraph_vector_scale(vector, -1.0);
      }
    }

    /* Correction for numeric inaccuracies (eliminating -0.0) */
    for (i=0; i<options->n; i++) {
      if (VECTOR(*vector)[i] < 0)
		    VECTOR(*vector)[i] = 0;
    }
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
 * \function igraph_eigenvector_centrality
 * Eigenvector centrality of the vertices
 * 
 * Eigenvector centrality is a measure of the importance of a node in a
 * network. It assigns relative scores to all nodes in the network based
 * on the principle that connections to high-scoring nodes contribute
 * more to the score of the node in question than equal connections to
 * low-scoring nodes. In practice, this is determined by calculating the
 * eigenvector corresponding to the largest positive eigenvalue of the
 * adjacency matrix. The centrality scores returned by igraph are always
 * normalized such that the largest eigenvector centrality score is one
 * (with one exception, see below).
 *
 * </para><para>
 * Since the eigenvector centrality scores of nodes in different components
 * do not affect each other, it may be beneficial for large graphs to
 * decompose it first into weakly connected components and calculate the
 * centrality scores individually for each component.
 *
 * </para><para>
 * Also note that the adjacency matrix of a directed acyclic graph or the
 * adjacency matrix of an empty graph does not possess positive eigenvalues,
 * therefore the eigenvector centrality is not defined for these graphs.
 * igraph will return an eigenvalue of zero in such cases. The eigenvector
 * centralities will all be equal for an empty graph and will all be zeros
 * for a directed acyclic graph. Such pathological cases can be detected
 * by asking igraph to calculate the eigenvalue as well (using the \p value
 * parameter, see below) and checking whether the eigenvalue is very close
 * to zero.
 *
 * \param graph The input graph. It might be directed.
 * \param vector Pointer to an initialized vector, it will be resized
 *     as needed. The result of the computation is stored here. It can
 *     be a null pointer, then it is ignored.
 * \param value If not a null pointer, then the eigenvalue
 *     corresponding to the found eigenvector is stored here.
 * \param directed Boolean scalar, whether to consider edge directions
 *     in a directed graph. It is ignored for undirected graphs.
 * \param scale If not zero then the result will be scaled such that
 *     the absolute value of the maximum centrality is one.
 * \param weights A null pointer (=no edge weights), or a vector
 *     giving the weights of the edges. The algorithm might result
 *     complex numbers is some weights are negative. In this case only
 *     the real part is reported.
 * \param options Options to ARPACK. See \ref igraph_arpack_options_t
 *    for details. Note that the function overwrites the
 *    <code>n</code> (number of vertices) parameter and 
 *    it always starts the calculation from a non-random vector
 *    calculated based on the degree of the vertices.
 * \return Error code.
 * 
 * Time complexity: depends on the input graph, usually it is O(|V|+|E|).
 * 
 * \sa \ref igraph_pagerank and \ref igraph_personalized_pagerank for 
 *   modifications of eigenvector centrality.
 * 
 * \example examples/simple/eigenvector_centrality.c
 */

int igraph_eigenvector_centrality(const igraph_t *graph, 
				  igraph_vector_t *vector,
				  igraph_real_t *value, 
				  igraph_bool_t directed, igraph_bool_t scale,
				  const igraph_vector_t *weights,
				  igraph_arpack_options_t *options) {

  if (directed && igraph_is_directed(graph)) {
    return igraph_eigenvector_centrality_directed(graph, vector, value,
						  scale, weights, options);
  } else {
    return igraph_eigenvector_centrality_undirected(graph, vector, value, 
						    scale, weights, options);
  }
}

/* struct for the unweighted variant of the HITS algorithm */
typedef struct igraph_i_kleinberg_data_t {
  igraph_adjlist_t *in;
  igraph_adjlist_t *out;
  igraph_vector_t *tmp;
} igraph_i_kleinberg_data_t;

/* struct for the weighted variant of the HITS algorithm */
typedef struct igraph_i_kleinberg_data2_t {
  const igraph_t *graph;
  igraph_inclist_t *in;
  igraph_inclist_t *out;
  igraph_vector_t *tmp;
  const igraph_vector_t *weights;
} igraph_i_kleinberg_data2_t;

/* ARPACK auxiliary routine for the unweighted HITS algorithm */
int igraph_i_kleinberg_unweighted(igraph_real_t *to,
                                  const igraph_real_t *from,
                                  int n, void *extra) {
  igraph_i_kleinberg_data_t *data = (igraph_i_kleinberg_data_t*)extra;
  igraph_adjlist_t *in = data->in;
  igraph_adjlist_t *out = data->out;
  igraph_vector_t *tmp = data->tmp;
  igraph_vector_int_t *neis;
  long int i, j, nlen;
  
  for (i=0; i<n; i++) {
    neis=igraph_adjlist_get(in, i);
    nlen=igraph_vector_int_size(neis);
    VECTOR(*tmp)[i]=0.0;
    for (j=0; j<nlen; j++) {
      long int nei=(long int) VECTOR(*neis)[j];
      VECTOR(*tmp)[i] += from[nei];
    }
  }
  
  for (i=0; i<n; i++) {
    neis=igraph_adjlist_get(out, i);
    nlen=igraph_vector_int_size(neis);
    to[i]=0.0;
    for (j=0; j<nlen; j++) {
      long int nei=(long int) VECTOR(*neis)[j];
      to[i] += VECTOR(*tmp)[nei];
    }
  }      
  
  return 0;
}

/* ARPACK auxiliary routine for the weighted HITS algorithm */
int igraph_i_kleinberg_weighted(igraph_real_t *to,
                                const igraph_real_t *from,
                                int n, void *extra) {

  igraph_i_kleinberg_data2_t *data = (igraph_i_kleinberg_data2_t*)extra;
  igraph_inclist_t *in = data->in; 
  igraph_inclist_t *out = data->out; 
  igraph_vector_t *tmp = data->tmp;
  const igraph_vector_t *weights = data->weights; 
  const igraph_t *g = data->graph;
  igraph_vector_int_t *neis;
  long int i, j, nlen;
  
  for (i=0; i<n; i++) {
    neis=igraph_inclist_get(in, i);
    nlen=igraph_vector_int_size(neis);
    VECTOR(*tmp)[i]=0.0;
    for (j=0; j<nlen; j++) {
      long int nei_edge = (long int) VECTOR(*neis)[j];
      long int nei=IGRAPH_OTHER(g, nei_edge, i);
      VECTOR(*tmp)[i] += from[nei] * VECTOR(*weights)[nei_edge];
    }
  }
  
  for (i=0; i<n; i++) {
    neis=igraph_inclist_get(out, i);
    nlen=igraph_vector_int_size(neis);
    to[i]=0.0;
    for (j=0; j<nlen; j++) {
      long int nei_edge=(long int) VECTOR(*neis)[j];
      long int nei=IGRAPH_OTHER(g, nei_edge, i);
      to[i] += VECTOR(*tmp)[nei] * VECTOR(*weights)[nei_edge];
    }
  }      
  
  return 0;
}

int igraph_i_kleinberg(const igraph_t *graph, igraph_vector_t *vector,
		       igraph_real_t *value, igraph_bool_t scale,
			   const igraph_vector_t *weights,
		       igraph_arpack_options_t *options, int inout) {
  
  igraph_adjlist_t myinadjlist, myoutadjlist;
  igraph_inclist_t myininclist, myoutinclist;
  igraph_adjlist_t *inadjlist, *outadjlist;
  igraph_inclist_t *ininclist, *outinclist;
  igraph_vector_t tmp;
  igraph_vector_t values;
  igraph_matrix_t vectors;
  igraph_i_kleinberg_data_t extra;
  igraph_i_kleinberg_data2_t extra2;
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

  if (weights) {
	igraph_real_t min, max;

    if (igraph_vector_size(weights) != igraph_ecount(graph)) {
      IGRAPH_ERROR("Invalid length of weights vector when calculating "
                   "hub or authority scores", IGRAPH_EINVAL);
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

  options->n=igraph_vcount(graph);
  options->start=1;     /* no random start vector */
  
  IGRAPH_VECTOR_INIT_FINALLY(&values, 0);
  IGRAPH_MATRIX_INIT_FINALLY(&vectors, options->n, 1);
  IGRAPH_VECTOR_INIT_FINALLY(&tmp, options->n);
  
  if (inout==0) {
    inadjlist=&myinadjlist; 
    outadjlist=&myoutadjlist;
    ininclist=&myininclist;
    outinclist=&myoutinclist;
  } else if (inout==1) {
    inadjlist=&myoutadjlist;
    outadjlist=&myinadjlist;
    ininclist=&myoutinclist;
    outinclist=&myininclist;
  } else {
    /* This should not happen */
    IGRAPH_ERROR("Invalid 'inout' argument, please do not call "
                 "this function directly", IGRAPH_FAILURE);
  }

  if (weights == 0) {
    IGRAPH_CHECK(igraph_adjlist_init(graph, &myinadjlist, IGRAPH_IN));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &myinadjlist);
    IGRAPH_CHECK(igraph_adjlist_init(graph, &myoutadjlist, IGRAPH_OUT));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &myoutadjlist);
  } else {
    IGRAPH_CHECK(igraph_inclist_init(graph, &myininclist, IGRAPH_IN));
    IGRAPH_FINALLY(igraph_inclist_destroy, &myininclist);
    IGRAPH_CHECK(igraph_inclist_init(graph, &myoutinclist, IGRAPH_OUT));
    IGRAPH_FINALLY(igraph_inclist_destroy, &myoutinclist);
  }

  IGRAPH_CHECK(igraph_degree(graph, &tmp, igraph_vss_all(), IGRAPH_ALL, 0));
  for (i=0; i<options->n; i++) {
    if (VECTOR(tmp)[i] != 0) { 
      MATRIX(vectors, i, 0) = VECTOR(tmp)[i];
    } else {
      MATRIX(vectors, i, 0) = 1.0;
    }
  }
	
  extra.in=inadjlist; extra.out=outadjlist; extra.tmp=&tmp;
  extra2.in=ininclist; extra2.out=outinclist; extra2.tmp=&tmp;
  extra2.graph=graph; extra2.weights=weights;

  options->nev = 1;
  options->ncv = 0;   /* 0 means "automatic" in igraph_arpack_rssolve */
  options->which[0]='L'; options->which[1]='M';

  if (weights == 0) {
    IGRAPH_CHECK(igraph_arpack_rssolve(igraph_i_kleinberg_unweighted, &extra,
                                       options, 0, &values, &vectors));
    igraph_adjlist_destroy(&myoutadjlist);
    igraph_adjlist_destroy(&myinadjlist);
	IGRAPH_FINALLY_CLEAN(2);
  } else {
    IGRAPH_CHECK(igraph_arpack_rssolve(igraph_i_kleinberg_weighted, &extra2,
                                       options, 0, &values, &vectors));
    igraph_inclist_destroy(&myoutinclist);
    igraph_inclist_destroy(&myininclist);
	IGRAPH_FINALLY_CLEAN(2);
  }

  igraph_vector_destroy(&tmp);
  IGRAPH_FINALLY_CLEAN(1);

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
    if (scale && amax!=0) {
      igraph_vector_scale(vector, 1/VECTOR(*vector)[which]);
    } else if (igraph_i_vector_mostly_negative(vector)) {
      igraph_vector_scale(vector, -1.0);
	}

    /* Correction for numeric inaccuracies (eliminating -0.0) */
    for (i=0; i<options->n; i++) {
      if (VECTOR(*vector)[i] < 0)
		  VECTOR(*vector)[i] = 0;
    }
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
 * \param scale If not zero then the result will be scaled such that
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
 * \sa \ref igraph_authority_score() for the companion measure,
 * \ref igraph_pagerank(), \ref igraph_personalized_pagerank(),
 * \ref igraph_eigenvector_centrality() for similar measures.
 */

int igraph_hub_score(const igraph_t *graph, igraph_vector_t *vector,
		     igraph_real_t *value, igraph_bool_t scale,
			 const igraph_vector_t *weights,
		     igraph_arpack_options_t *options) {

  return igraph_i_kleinberg(graph, vector, value, scale, weights, options, 0);
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
 * \param scale If not zero then the result will be scaled such that
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
 * \sa \ref igraph_hub_score() for the companion measure,
 * \ref igraph_pagerank(), \ref igraph_personalized_pagerank(),
 * \ref igraph_eigenvector_centrality() for similar measures.
 */
			    
int igraph_authority_score(const igraph_t *graph, igraph_vector_t *vector,
			   igraph_real_t *value, igraph_bool_t scale,
			   const igraph_vector_t *weights,
			   igraph_arpack_options_t *options) {

  return igraph_i_kleinberg(graph, vector, value, scale, weights, options, 1);
}

typedef struct igraph_i_pagerank_data_t {
  const igraph_t *graph;
  igraph_adjlist_t *adjlist;
  igraph_real_t damping;
  igraph_vector_t *outdegree;
  igraph_vector_t *tmp;
  igraph_vector_t *reset;
} igraph_i_pagerank_data_t;

typedef struct igraph_i_pagerank_data2_t {
  const igraph_t *graph;
  igraph_inclist_t *inclist;
  const igraph_vector_t *weights;
  igraph_real_t damping;
  igraph_vector_t *outdegree;
  igraph_vector_t *tmp;
  igraph_vector_t *reset;
} igraph_i_pagerank_data2_t;

int igraph_i_pagerank(igraph_real_t *to, const igraph_real_t *from,
		      int n, void *extra) {
  
  igraph_i_pagerank_data_t *data=extra;
  igraph_adjlist_t *adjlist=data->adjlist;
  igraph_vector_t *outdegree=data->outdegree;
  igraph_vector_t *tmp=data->tmp;
  igraph_vector_t *reset=data->reset;
  igraph_vector_int_t *neis;
  long int i, j, nlen;
  igraph_real_t sumfrom=0.0;
  igraph_real_t fact=1-data->damping;

  /* Calculate p(x) / outdegree(x) in advance for all the vertices.
   * Note that we may divide by zero here; this is intentional since
   * we won't use those values and we save a comparison this way.
   * At the same time, we calculate the global probability of a
   * random jump in `sumfrom`. For vertices with no outgoing edges,
   * we will surely jump from there if we are there, hence those
   * vertices contribute p(x) to the teleportation probability.
   * For vertices with some outgoing edges, we jump from there with
   * probability `fact` if we are there, hence they contribute
   * p(x)*fact */
  for (i=0; i<n; i++) {
    sumfrom += VECTOR(*outdegree)[i]!=0 ? from[i] * fact : from[i];
    VECTOR(*tmp)[i] = from[i] / VECTOR(*outdegree)[i];
  }

  /* Here we calculate the part of the `to` vector that results from
   * moving along links (and not from teleportation) */
  for (i=0; i<n; i++) {
    neis=igraph_adjlist_get(adjlist, i);
    nlen=igraph_vector_int_size(neis);
    to[i]=0.0;
    for (j=0; j<nlen; j++) {
      long int nei=(long int) VECTOR(*neis)[j];
      to[i] += VECTOR(*tmp)[nei];
    }
    to[i] *= data->damping;
  }

  /* Now we add the contribution from random jumps. `reset` is a vector
   * that defines the probability of ending up in vertex i after a jump.
   * `sumfrom` is the global probability of jumping as mentioned above. */
  /* printf("sumfrom = %.6f\n", (float)sumfrom); */

  if (reset) {
    /* Running personalized PageRank */
    for (i=0; i<n; i++) {
      to[i] += sumfrom * VECTOR(*reset)[i];
    }
  } else {
    /* Traditional PageRank with uniform reset vector */
    sumfrom /= n;
    for (i=0; i<n; i++) {
      to[i] += sumfrom;
    }
  }

  return 0;
}

int igraph_i_pagerank2(igraph_real_t *to, const igraph_real_t *from,
		       int n, void *extra) {

  igraph_i_pagerank_data2_t *data=extra;
  const igraph_t *graph=data->graph;
  igraph_inclist_t *inclist=data->inclist;
  const igraph_vector_t *weights=data->weights;
  igraph_vector_t *outdegree=data->outdegree;
  igraph_vector_t *tmp=data->tmp;
  igraph_vector_t *reset=data->reset;
  long int i, j, nlen;
  igraph_real_t sumfrom=0.0;
  igraph_vector_int_t *neis;
  igraph_real_t fact=1-data->damping;

  /*
  printf("PageRank weighted: multiplying vector: ");
  for (i=0; i<n; i++) { printf(" %.4f", from[i]); }
  printf("\n");
  */

  for (i=0; i<n; i++) {
    sumfrom += VECTOR(*outdegree)[i]!=0 ? from[i] * fact : from[i];
    VECTOR(*tmp)[i] = from[i] / VECTOR(*outdegree)[i];
  }
  
  for (i=0; i<n; i++) {
    neis=igraph_inclist_get(inclist, i);
    nlen=igraph_vector_int_size(neis);
    to[i]=0.0;
    for (j=0; j<nlen; j++) {
      long int edge=(long int) VECTOR(*neis)[j];
      long int nei=IGRAPH_OTHER(graph, edge, i);
      to[i] += VECTOR(*weights)[edge] * VECTOR(*tmp)[nei];
    }
    to[i] *= data->damping;
  }

  /* printf("sumfrom = %.6f\n", (float)sumfrom); */

  if (reset) {
    /* Running personalized PageRank */
    for (i=0; i<n; i++) {
      to[i] += sumfrom * VECTOR(*reset)[i];
    }
  } else {
    /* Traditional PageRank with uniform reset vector */
    sumfrom /= n;
    for (i=0; i<n; i++) {
      to[i] += sumfrom;
    }
  }

  /*
  printf("PageRank weighted: multiplied vector: ");
  for (i=0; i<n; i++) { printf(" %.4f", to[i]); }
  printf("\n");
  */

  return 0;
}

/**
 * \function igraph_pagerank
 * \brief Calculates the Google PageRank for the specified vertices.
 *
 * Starting from version 0.7, igraph has three PageRank implementations,
 * and the user can choose between them. The first implementation is
 * \c IGRAPH_PAGERANK_ALGO_POWER, also available as the (now
 * deprecated) function \ref igraph_pagerank_old(). The second
 * implementation is based on the ARPACK library, this was the default
 * before igraph version 0.7: \c IGRAPH_PAGERANK_ALGO_ARPACK.
 *
 * The third and recommmended implementation is \c
 * IGRAPH_PAGERANK_ALGO_PRPACK. This is using the the PRPACK package,
 * see https://github.com/dgleich/prpack .
 *
 * </para><para>
 * Please note that the PageRank of a given vertex depends on the PageRank
 * of all other vertices, so even if you want to calculate the PageRank for
 * only some of the vertices, all of them must be calculated. Requesting
 * the PageRank for only some of the vertices does not result in any
 * performance increase at all.
 * </para>
 * 
 * <para>
 * For the explanation of the PageRank algorithm, see the following
 * webpage:
 * http://infolab.stanford.edu/~backrub/google.html , or the
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
 * \param algo The PageRank implementation to use. Possible values:
 *    \c IGRAPH_PAGERANK_ALGO_POWER, \c IGRAPH_PAGERANK_ALGO_ARPACK,
 *    \c IGRAPH_PAGERANK_ALGO_PRPACK.
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
 * \param options Options to the power method or ARPACK. For the power
 *    method, \c IGRAPH_PAGERANK_ALGO_POWER it must be a pointer to
 *    a \ref igraph_pagerank_power_options_t object.
 *    For \c IGRAPH_PAGERANK_ALGO_ARPACK it must be a pointer to an
 *    \ref igraph_arpack_options_t object. See \ref igraph_arpack_options_t
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
 * Time complexity: depends on the input graph, usually it is O(|E|),
 * the number of edges.
 * 
 * \sa \ref igraph_pagerank_old() for the old implementation, 
 * \ref igraph_personalized_pagerank() and \ref igraph_personalized_pagerank_vs()
 * for the personalized PageRank measure, \ref igraph_arpack_rssolve() and
 * \ref igraph_arpack_rnsolve() for the underlying machinery.
 * 
 * \example examples/simple/igraph_pagerank.c
 */

int igraph_pagerank(const igraph_t *graph, igraph_pagerank_algo_t algo,
		    igraph_vector_t *vector,
		    igraph_real_t *value, const igraph_vs_t vids,
		    igraph_bool_t directed, igraph_real_t damping, 
		    const igraph_vector_t *weights, void *options) {
  return igraph_personalized_pagerank(graph, algo, vector, value, vids, 
				      directed, damping, 0, weights, 
				      options);
}

/**
 * \function igraph_personalized_pagerank_vs
 * \brief Calculates the personalized Google PageRank for the specified vertices.
 * 
 * The personalized PageRank is similar to the original PageRank measure, but the
 * random walk is reset in every step with probability 1-damping to a non-uniform
 * distribution (instead of the uniform distribution in the original PageRank measure.
 *
 * </para><para>
 * This simplified interface takes a vertex sequence and resets the random walk to
 * one of the vertices in the specified vertex sequence, chosen uniformly. A typical
 * application of personalized PageRank is when the random walk is reset to the same
 * vertex every time - this can easily be achieved using \ref igraph_vss_1() which
 * generates a vertex sequence containing only a single vertex.
 * 
 * </para><para>
 * Please note that the personalized PageRank of a given vertex depends on the
 * personalized PageRank of all other vertices, so even if you want to calculate
 * the personalized PageRank for only some of the vertices, all of them must be
 * calculated. Requesting the personalized PageRank for only some of the vertices
 * does not result in any performance increase at all.
 * </para>
 * 
 * <para>
 * \param graph The graph object.
 * \param algo The PageRank implementation to use. Possible values:
 *    \c IGRAPH_PAGERANK_ALGO_POWER, \c IGRAPH_PAGERANK_ALGO_ARPACK,
 *    \c IGRAPH_PAGERANK_ALGO_PRPACK.
 * \param vector Pointer to an initialized vector, the result is
 *    stored here. It is resized as needed.
 * \param value Pointer to a real variable, the eigenvalue
 *    corresponding to the PageRank vector is stored here. It should
 *    be always exactly one.
 * \param vids The vertex ids for which the PageRank is returned.
 * \param directed Boolean, whether to consider the directedness of
 *    the edges. This is ignored for undirected graphs.
 * \param damping The damping factor ("d" in the original paper)
 * \param reset_vids IDs of the vertices used when resetting the random walk.
 * \param weights Optional edge weights, it is either a null pointer,
 *    then the edges are not weighted, or a vector of the same length
 *    as the number of edges.
 * \param options Options to the power method or ARPACK. For the power
 *    method, \c IGRAPH_PAGERANK_ALGO_POWER it must be a pointer to
 *    a \ref igraph_pagerank_power_options_t object.
 *    For \c IGRAPH_PAGERANK_ALGO_ARPACK it must be a pointer to an
 *    \ref igraph_arpack_options_t object. See \ref igraph_arpack_options_t
 *    for details. Note that the function overwrites the
 *    <code>n</code> (number of vertices), <code>nev</code> (1),
 *    <code>ncv</code> (3) and <code>which</code> (LM) parameters and
 *    it always starts the calculation from a non-random vector
 *    calculated based on the degree of the vertices.
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory for
 *         temporary data. 
 *         \c IGRAPH_EINVVID, invalid vertex id in
 *         \p vids or an empty reset vertex sequence in
 *         \p vids_reset.
 * 
 * Time complexity: depends on the input graph, usually it is O(|E|),
 * the number of edges.
 * 
 * \sa \ref igraph_pagerank() for the non-personalized implementation, 
 * \ref igraph_arpack_rssolve() and \ref igraph_arpack_rnsolve() for
 * the underlying machinery.
 */

int igraph_personalized_pagerank_vs(const igraph_t *graph, 
		    igraph_pagerank_algo_t algo, igraph_vector_t *vector,
		    igraph_real_t *value, const igraph_vs_t vids,
		    igraph_bool_t directed, igraph_real_t damping, 
		    igraph_vs_t reset_vids,
		    const igraph_vector_t *weights,
		    void *options) {
	igraph_vector_t reset;
	igraph_vit_t vit;

	IGRAPH_VECTOR_INIT_FINALLY(&reset, igraph_vcount(graph));
	IGRAPH_CHECK(igraph_vit_create(graph, reset_vids, &vit));
	IGRAPH_FINALLY(igraph_vit_destroy, &vit);

	while (!IGRAPH_VIT_END(vit)) {
		VECTOR(reset)[(long int)IGRAPH_VIT_GET(vit)]++;
		IGRAPH_VIT_NEXT(vit);
	}
	igraph_vit_destroy(&vit);
	IGRAPH_FINALLY_CLEAN(1);
	
	IGRAPH_CHECK(igraph_personalized_pagerank(graph, algo, vector, 
						  value, vids, directed, 
						  damping, &reset, weights,
						  options));

	igraph_vector_destroy(&reset);
	IGRAPH_FINALLY_CLEAN(1);

	return 0;
}

/**
 * \function igraph_personalized_pagerank
 * \brief Calculates the personalized Google PageRank for the specified vertices.
 * 
 * The personalized PageRank is similar to the original PageRank measure, but the
 * random walk is reset in every step with probability 1-damping to a non-uniform
 * distribution (instead of the uniform distribution in the original PageRank measure.
 *
 * </para><para>
 * Please note that the personalized PageRank of a given vertex depends on the
 * personalized PageRank of all other vertices, so even if you want to calculate
 * the personalized PageRank for only some of the vertices, all of them must be
 * calculated. Requesting the personalized PageRank for only some of the vertices
 * does not result in any performance increase at all.
 * </para>
 * 
 * <para>
 * \param graph The graph object.
 * \param algo The PageRank implementation to use. Possible values:
 *    \c IGRAPH_PAGERANK_ALGO_POWER, \c IGRAPH_PAGERANK_ALGO_ARPACK,
 *    \c IGRAPH_PAGERANK_ALGO_PRPACK.
 * \param vector Pointer to an initialized vector, the result is
 *    stored here. It is resized as needed.
 * \param value Pointer to a real variable, the eigenvalue
 *    corresponding to the PageRank vector is stored here. It should
 *    be always exactly one.
 * \param vids The vertex ids for which the PageRank is returned.
 * \param directed Boolean, whether to consider the directedness of
 *    the edges. This is ignored for undirected graphs.
 * \param damping The damping factor ("d" in the original paper)
 * \param reset The probability distribution over the vertices used when
 *    resetting the random walk. It is either a null pointer (denoting
 *    a uniform choice that results in the original PageRank measure)
 *    or a vector of the same length as the number of vertices.
 * \param weights Optional edge weights, it is either a null pointer,
 *    then the edges are not weighted, or a vector of the same length
 *    as the number of edges.
 * \param options Options to the power method or ARPACK. For the power
 *    method, \c IGRAPH_PAGERANK_ALGO_POWER it must be a pointer to
 *    a \ref igraph_pagerank_power_options_t object.
 *    For \c IGRAPH_PAGERANK_ALGO_ARPACK it must be a pointer to an
 *    \ref igraph_arpack_options_t object. See \ref igraph_arpack_options_t
 *    for details. Note that the function overwrites the
 *    <code>n</code> (number of vertices), <code>nev</code> (1),
 *    <code>ncv</code> (3) and <code>which</code> (LM) parameters and
 *    it always starts the calculation from a non-random vector
 *    calculated based on the degree of the vertices.
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory for
 *         temporary data. 
 *         \c IGRAPH_EINVVID, invalid vertex id in
 *         \p vids or an invalid reset vector in \p reset. 
 * 
 * Time complexity: depends on the input graph, usually it is O(|E|),
 * the number of edges.
 * 
 * \sa \ref igraph_pagerank() for the non-personalized implementation, 
 * \ref igraph_arpack_rssolve() and \ref igraph_arpack_rnsolve() for
 * the underlying machinery.
 */
int igraph_personalized_pagerank(const igraph_t *graph, 
		    igraph_pagerank_algo_t algo, igraph_vector_t *vector,
		    igraph_real_t *value, const igraph_vs_t vids,
		    igraph_bool_t directed, igraph_real_t damping, 
		    igraph_vector_t *reset,
		    const igraph_vector_t *weights,
		    void *options) {

  if (algo == IGRAPH_PAGERANK_ALGO_POWER) {
    igraph_pagerank_power_options_t *o = 
      (igraph_pagerank_power_options_t *) options;
    if (reset) { 
      IGRAPH_WARNING("Cannot use weights with power method, "
		     "weights will be ignored");
    }
    return igraph_pagerank_old(graph, vector, vids, directed, 
			       o->niter, o->eps, damping, 
			       /*old=*/ 0);
  } else if (algo == IGRAPH_PAGERANK_ALGO_ARPACK) {
    igraph_arpack_options_t *o= (igraph_arpack_options_t*) options;
    return igraph_personalized_pagerank_arpack(graph, vector, value, vids,
					       directed, damping, reset, 
					       weights, o);
  } else if (algo == IGRAPH_PAGERANK_ALGO_PRPACK) {
    return igraph_personalized_pagerank_prpack(graph, vector, value, vids,
					       directed, damping, reset, 
					       weights);
  } else {
    IGRAPH_ERROR("Unknown PageRank algorithm", IGRAPH_EINVAL);
  }

  return 0;
}

/*
 * ARPACK-based implementation of \c igraph_personalized_pagerank.
 *
 * See \c igraph_personalized_pagerank for the documentation of the parameters.
 */
int igraph_personalized_pagerank_arpack(const igraph_t *graph, igraph_vector_t *vector,
		    igraph_real_t *value, const igraph_vs_t vids,
		    igraph_bool_t directed, igraph_real_t damping, 
		    igraph_vector_t *reset,
		    const igraph_vector_t *weights,
		    igraph_arpack_options_t *options) {
  igraph_matrix_t values;
  igraph_matrix_t vectors;
  igraph_neimode_t dirmode;
  igraph_vector_t outdegree;
  igraph_vector_t indegree;
  igraph_vector_t tmp;

  long int i;
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);

  if (no_of_edges == 0) {
    /* special case: empty graph */
    if (value)
      *value = 1.0;
    if (vector) {
      igraph_vector_resize(vector, no_of_nodes);
      igraph_vector_fill(vector, 1.0 / no_of_nodes);
    }
	return IGRAPH_SUCCESS;
  }

  options->n = (int) no_of_nodes;
  options->nev = 1;
  options->ncv = 0;   /* 0 means "automatic" in igraph_arpack_rnsolve */
  options->which[0]='L'; options->which[1]='M';
  options->start = 1;		/* no random start vector */

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
        *value = 1.0;
      if (vector) {
        igraph_vector_resize(vector, igraph_vcount(graph));
        igraph_vector_fill(vector, 1.0 / no_of_nodes);
      }
      return IGRAPH_SUCCESS;
    }
  }

  if (reset && igraph_vector_size(reset) != no_of_nodes)
  {
    IGRAPH_ERROR("Invalid length of reset vector when calculating "
		 "personalized PageRank scores", IGRAPH_EINVAL);
  }

  IGRAPH_MATRIX_INIT_FINALLY(&values, 0, 0);
  IGRAPH_MATRIX_INIT_FINALLY(&vectors, options->n, 1);

  if (directed) { dirmode=IGRAPH_IN; } else { dirmode=IGRAPH_ALL; }

  IGRAPH_VECTOR_INIT_FINALLY(&indegree, options->n);
  IGRAPH_VECTOR_INIT_FINALLY(&outdegree, options->n);
  IGRAPH_VECTOR_INIT_FINALLY(&tmp, options->n);

  RNG_BEGIN();

  if (reset) {
	/* Normalize reset vector so the sum is 1 */
	double reset_sum;
	if (igraph_vector_min(reset) < 0)
	  IGRAPH_ERROR("the reset vector must not contain negative elements", IGRAPH_EINVAL);
	reset_sum = igraph_vector_sum(reset);
	if (reset_sum == 0)
	  IGRAPH_ERROR("the sum of the elements in the reset vector must not be zero", IGRAPH_EINVAL);
	igraph_vector_scale(reset, 1.0/reset_sum);
  }

  if (!weights) {
    
    igraph_adjlist_t adjlist;
    igraph_i_pagerank_data_t data = { graph, &adjlist, damping,
				      &outdegree, &tmp, reset };

    IGRAPH_CHECK(igraph_degree(graph, &outdegree, igraph_vss_all(),
			       directed ? IGRAPH_OUT : IGRAPH_ALL, /*loops=*/ 0));
    IGRAPH_CHECK(igraph_degree(graph, &indegree, igraph_vss_all(),
			       directed ? IGRAPH_IN : IGRAPH_ALL, /*loops=*/ 0));
    /* Set up an appropriate starting vector. We start from the in-degrees
     * plus some small random noise to avoid convergence problems */
    for (i=0; i<options->n; i++) {
      if (VECTOR(indegree)[i])
        MATRIX(vectors, i, 0) = VECTOR(indegree)[i] + RNG_UNIF(-1e-4, 1e-4);
      else
        MATRIX(vectors, i, 0) = 1;
    }

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, dirmode));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    IGRAPH_CHECK(igraph_arpack_rnsolve(igraph_i_pagerank,
				       &data, options, 0, &values, &vectors));

    igraph_adjlist_destroy(&adjlist);
    IGRAPH_FINALLY_CLEAN(1);
    
  } else {
    
    igraph_inclist_t inclist;
    igraph_bool_t negative_weight_warned = 0;
    igraph_i_pagerank_data2_t data = { graph, &inclist, weights,
				       damping, &outdegree, &tmp, reset };    

    IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, dirmode));
    IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);

    /* Weighted degree */
    for (i=0; i<no_of_edges; i++) {
      long int from=IGRAPH_FROM(graph, i);
      long int to=IGRAPH_TO(graph, i);
      igraph_real_t weight=VECTOR(*weights)[i];
      if (weight < 0 && !negative_weight_warned) {
        IGRAPH_WARNING("replacing negative weights with zeros");
        weight = 0;
        negative_weight_warned = 1;
      }
      VECTOR(outdegree)[from] += weight;
      VECTOR(indegree) [to]   += weight;
      if (!directed) { 
        VECTOR(outdegree)[to]   += weight;
        VECTOR(indegree) [from] += weight;
      }
    }
    /* Set up an appropriate starting vector. We start from the in-degrees
     * plus some small random noise to avoid convergence problems */
    for (i=0; i<options->n; i++) {
      if (VECTOR(indegree)[i])
        MATRIX(vectors, i, 0) = VECTOR(indegree)[i] + RNG_UNIF(-1e-4, 1e-4);
      else
        MATRIX(vectors, i, 0) = 1;
    }
    
    IGRAPH_CHECK(igraph_arpack_rnsolve(igraph_i_pagerank2,
				       &data, options, 0, &values, &vectors));
    
    igraph_inclist_destroy(&inclist);
    IGRAPH_FINALLY_CLEAN(1);
  }

  RNG_END();

  igraph_vector_destroy(&tmp);
  igraph_vector_destroy(&outdegree);
  igraph_vector_destroy(&indegree);
  IGRAPH_FINALLY_CLEAN(3);

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
 * \param weights An optional vector containing edge weights for 
 *        calculating weighted betweenness. Supply a null pointer here
 *        for unweighted betweenness.
 * \param nobigint Logical, if true, then we don't use big integers
 *        for the calculation, setting this to 1 (=true) should
 *        work for most graphs. It is currently ignored for weighted
 *        graphs.
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
 * 
 * \example examples/simple/igraph_betweenness.c
 */
int igraph_betweenness(const igraph_t *graph, igraph_vector_t *res,
		       const igraph_vs_t vids, igraph_bool_t directed, 
		       const igraph_vector_t* weights, igraph_bool_t nobigint) {
  return igraph_betweenness_estimate(graph, res, vids, directed, -1, weights,
				     nobigint);
}

int igraph_i_betweenness_estimate_weighted(const igraph_t *graph, 
					 igraph_vector_t *res, 
					 const igraph_vs_t vids, 
					 igraph_bool_t directed,
					 igraph_real_t cutoff, 
					 const igraph_vector_t *weights, 
					 igraph_bool_t nobigint) {

  igraph_real_t minweight;
  igraph_integer_t no_of_nodes=(igraph_integer_t) igraph_vcount(graph);
  igraph_integer_t no_of_edges=(igraph_integer_t) igraph_ecount(graph);
  igraph_2wheap_t Q;
  igraph_inclist_t inclist;
  igraph_adjlist_t fathers;
  long int source, j;
  igraph_stack_t S;
  igraph_neimode_t mode= directed ? IGRAPH_OUT : IGRAPH_ALL;
  igraph_vector_t dist, nrgeo, tmpscore;
  igraph_vector_t v_tmpres, *tmpres=&v_tmpres;
  igraph_vit_t vit;
  int cmp_result;
  const double eps = IGRAPH_SHORTEST_PATH_EPSILON;

  IGRAPH_UNUSED(nobigint);

  if (igraph_vector_size(weights) != no_of_edges) {
    IGRAPH_ERROR("Weight vector length does not match", IGRAPH_EINVAL);
  }
  minweight = igraph_vector_min(weights);
  if (minweight <= 0) {
    IGRAPH_ERROR("Weight vector must be positive", IGRAPH_EINVAL);
  }
  else if (minweight <= eps) {
    IGRAPH_WARNING("Some weights are smaller than epsilon, calculations may suffer from numerical precision.");
  }

  IGRAPH_CHECK(igraph_2wheap_init(&Q, no_of_nodes));
  IGRAPH_FINALLY(igraph_2wheap_destroy, &Q);
  IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, mode));  
  IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
  IGRAPH_CHECK(igraph_adjlist_init_empty(&fathers, no_of_nodes));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &fathers);

  IGRAPH_CHECK(igraph_stack_init(&S, no_of_nodes));
  IGRAPH_FINALLY(igraph_stack_destroy, &S);
  IGRAPH_VECTOR_INIT_FINALLY(&dist, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&tmpscore, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&nrgeo, no_of_nodes);

  if (igraph_vs_is_all(&vids)) {
    IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
    igraph_vector_null(res);
    tmpres=res;
  } else {
    IGRAPH_VECTOR_INIT_FINALLY(tmpres, no_of_nodes);
  }

  for (source=0; source<no_of_nodes; source++) {
    IGRAPH_PROGRESS("Betweenness centrality: ", 100.0*source/no_of_nodes, 0);
    IGRAPH_ALLOW_INTERRUPTION();

    igraph_2wheap_push_with_index(&Q, source, -1.0);
    VECTOR(dist)[source]=1.0;
    VECTOR(nrgeo)[source]=1;
    
    while (!igraph_2wheap_empty(&Q)) {
      long int minnei=igraph_2wheap_max_index(&Q);
      igraph_real_t mindist=-igraph_2wheap_delete_max(&Q);
      igraph_vector_int_t *neis;
      long int nlen;
      
      igraph_stack_push(&S, minnei);
      if (cutoff > 0 && VECTOR(dist)[minnei] >= cutoff+1.0) { continue; }
      
      /* Now check all neighbors of 'minnei' for a shorter path */
      neis=igraph_inclist_get(&inclist, minnei);
      nlen=igraph_vector_int_size(neis);
      for (j=0; j<nlen; j++) {
	long int edge=(long int) VECTOR(*neis)[j];
	long int to=IGRAPH_OTHER(graph, edge, minnei);
	igraph_real_t altdist=mindist + VECTOR(*weights)[edge];
	igraph_real_t curdist=VECTOR(dist)[to];

	if (curdist == 0) {
	  /* this means curdist is infinity */
	  cmp_result = -1;
	} else {
	  cmp_result = igraph_cmp_epsilon(altdist, curdist, eps);
	}
	
	if (curdist==0) {
	  /* This is the first non-infinite distance */
	  igraph_vector_int_t *v=igraph_adjlist_get(&fathers, to);
	  igraph_vector_int_resize(v,1);
	  VECTOR(*v)[0]=minnei;
	  VECTOR(nrgeo)[to] = VECTOR(nrgeo)[minnei];

	  VECTOR(dist)[to]=altdist;
	  IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, to, -altdist));
	} else if (cmp_result < 0) {
	  /* This is a shorter path */
	  igraph_vector_int_t *v=igraph_adjlist_get(&fathers, to);
	  igraph_vector_int_resize(v,1);
	  VECTOR(*v)[0]=minnei;
	  VECTOR(nrgeo)[to] = VECTOR(nrgeo)[minnei];

	  VECTOR(dist)[to]=altdist;
	  IGRAPH_CHECK(igraph_2wheap_modify(&Q, to, -altdist));
	} else if (cmp_result == 0) {
	  igraph_vector_int_t *v=igraph_adjlist_get(&fathers, to);
	  igraph_vector_int_push_back(v, minnei);
	  VECTOR(nrgeo)[to] += VECTOR(nrgeo)[minnei];
	}
      }
      
    } /* !igraph_2wheap_empty(&Q) */

    while (!igraph_stack_empty(&S)) {
      long int w=(long int) igraph_stack_pop(&S);
      igraph_vector_int_t *fatv=igraph_adjlist_get(&fathers, w);
      long int fatv_len=igraph_vector_int_size(fatv);
      for (j=0; j<fatv_len; j++) {
	long int f=(long int) VECTOR(*fatv)[j];
	VECTOR(tmpscore)[f] += VECTOR(nrgeo)[f]/VECTOR(nrgeo)[w] * (1+VECTOR(tmpscore)[w]);
      }
      if (w!=source) { VECTOR(*tmpres)[w] += VECTOR(tmpscore)[w]; }

      VECTOR(tmpscore)[w]=0;
      VECTOR(dist)[w]=0;
      VECTOR(nrgeo)[w]=0;
      igraph_vector_int_clear(igraph_adjlist_get(&fathers, w));
    }
    
  } /* source < no_of_nodes */

  if (!igraph_vs_is_all(&vids)) {
    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_VIT_SIZE(vit)));
    
    for (j=0, IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit);
	 IGRAPH_VIT_NEXT(vit), j++) {
      long int node=IGRAPH_VIT_GET(vit);
      VECTOR(*res)[j] = VECTOR(*tmpres)[node];
    }
    
    no_of_nodes = (igraph_integer_t) j;
    
    igraph_vit_destroy(&vit);
    igraph_vector_destroy(tmpres);
    IGRAPH_FINALLY_CLEAN(2);
  }

  if (!directed || !igraph_is_directed(graph)) {
    for (j=0; j<no_of_nodes; j++) {
      VECTOR(*res)[j] /= 2.0;
    }
  }
  
  IGRAPH_PROGRESS("Betweenness centrality: ", 100.0, 0);

  igraph_vector_destroy(&nrgeo);
  igraph_vector_destroy(&tmpscore);
  igraph_vector_destroy(&dist);
  igraph_stack_destroy(&S);
  igraph_adjlist_destroy(&fathers);
  igraph_inclist_destroy(&inclist);
  igraph_2wheap_destroy(&Q);
  IGRAPH_FINALLY_CLEAN(7);
  
  return 0;
}

void igraph_i_destroy_biguints(igraph_biguint_t *p) {
  igraph_biguint_t *p2 = p;
  while ( *((long int*)(p)) ) {
    igraph_biguint_destroy(p);
    p++;
  }
  igraph_Free(p2);
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
 * \param weights An optional vector containing edge weights for 
 *        calculating weighted betweenness. Supply a null pointer here
 *        for unweighted betweenness.
 * \param nobigint Logical, if true, then we don't use big integers
 *        for the calculation, setting this to 1 (=true) should
 *        work for most graphs. It is currently ignored for weighted
 *        graphs.
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
				igraph_real_t cutoff, 
				const igraph_vector_t *weights, 
				igraph_bool_t nobigint) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  long int *distance;
  unsigned long long int *nrgeo=0;  /* must be long long; consider grid
				       graphs for example */
  igraph_biguint_t *big_nrgeo=0;
  double *tmpscore;
  igraph_stack_t stack=IGRAPH_STACK_NULL;
  long int source;
  long int j, k, nneis;
  igraph_vector_int_t *neis;
  igraph_vector_t v_tmpres, *tmpres=&v_tmpres;
  igraph_vit_t vit;

  igraph_adjlist_t adjlist_out, adjlist_in;
  igraph_adjlist_t *adjlist_out_p, *adjlist_in_p;

  igraph_biguint_t D, R, T;

  if (weights) { 
    return igraph_i_betweenness_estimate_weighted(graph, res, vids, directed,
						cutoff, weights, nobigint);
  }

  if (!igraph_vs_is_all(&vids)) {
    /* subset */
    IGRAPH_VECTOR_INIT_FINALLY(tmpres, no_of_nodes);
  } else {
    /* only  */
    IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
    igraph_vector_null(res);
    tmpres=res;
  }

  directed=directed && igraph_is_directed(graph);
  if (directed) {
    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist_out, IGRAPH_OUT));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist_out);
    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist_in, IGRAPH_IN));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist_in);
    adjlist_out_p=&adjlist_out;
    adjlist_in_p=&adjlist_in;
  } else {
    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist_out, IGRAPH_ALL));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist_out);
    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist_in, IGRAPH_ALL));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist_in);
    adjlist_out_p=&adjlist_out;
    adjlist_in_p=&adjlist_in;
  }
  for (j=0; j<no_of_nodes; j++) {
    igraph_vector_int_clear(igraph_adjlist_get(adjlist_in_p, j));
  }
  
  distance=igraph_Calloc(no_of_nodes, long int);
  if (distance==0) {
    IGRAPH_ERROR("betweenness failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, distance);
  if (nobigint) {
    nrgeo=igraph_Calloc(no_of_nodes, unsigned long long int);
    if (nrgeo==0) {
      IGRAPH_ERROR("betweenness failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, nrgeo);
  } else {
    /* +1 is to have one containing zeros, when we free it, we stop
       at the zero */
    big_nrgeo=igraph_Calloc(no_of_nodes+1, igraph_biguint_t);
    if (!big_nrgeo) {
      IGRAPH_ERROR("betweenness failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_i_destroy_biguints, big_nrgeo);
    for (j=0; j<no_of_nodes; j++) {
      IGRAPH_CHECK(igraph_biguint_init(&big_nrgeo[j]));
    }
    IGRAPH_CHECK(igraph_biguint_init(&D));
    IGRAPH_FINALLY(igraph_biguint_destroy, &D);
    IGRAPH_CHECK(igraph_biguint_init(&R));
    IGRAPH_FINALLY(igraph_biguint_destroy, &R);
    IGRAPH_CHECK(igraph_biguint_init(&T));
    IGRAPH_FINALLY(igraph_biguint_destroy, &T);
  }
  tmpscore=igraph_Calloc(no_of_nodes, double);
  if (tmpscore==0) {
    IGRAPH_ERROR("betweenness failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, tmpscore);

  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
  igraph_stack_init(&stack, no_of_nodes);
  IGRAPH_FINALLY(igraph_stack_destroy, &stack);
    
  /* here we go */
  
  for (source=0; source<no_of_nodes; source++) {
    IGRAPH_PROGRESS("Betweenness centrality: ", 100.0*source/no_of_nodes, 0);
    IGRAPH_ALLOW_INTERRUPTION();

    IGRAPH_CHECK(igraph_dqueue_push(&q, source));
    if (nobigint) { 
      nrgeo[source]=1;
    } else {
      igraph_biguint_set_limb(&big_nrgeo[source], 1);
    }
    distance[source]=1;
    
    while (!igraph_dqueue_empty(&q)) {
      long int actnode=(long int) igraph_dqueue_pop(&q);
      IGRAPH_CHECK(igraph_stack_push(&stack, actnode));

      if (cutoff > 0 && distance[actnode] >= cutoff+1) { continue; }
      
      neis = igraph_adjlist_get(adjlist_out_p, actnode);
      nneis = igraph_vector_int_size(neis);
      for (j=0; j<nneis; j++) {
        long int neighbor=(long int) VECTOR(*neis)[j];
        if (distance[neighbor]==0) {
	  distance[neighbor]=distance[actnode]+1;
	  IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
	} 
	if (distance[neighbor]==distance[actnode]+1) {
	  igraph_vector_int_t *v=igraph_adjlist_get(adjlist_in_p, 
						    neighbor);
	  igraph_vector_int_push_back(v, actnode);
	  if (nobigint) { 
	    nrgeo[neighbor]+=nrgeo[actnode];
	  } else {
	    IGRAPH_CHECK(igraph_biguint_add(&big_nrgeo[neighbor],
					    &big_nrgeo[neighbor], 
					    &big_nrgeo[actnode]));
	  }
	}
      }
    } /* while !igraph_dqueue_empty */
    
    /* Ok, we've the distance of each node and also the number of
       shortest paths to them. Now we do an inverse search, starting
       with the farthest nodes. */
    while (!igraph_stack_empty(&stack)) {
      long int actnode=(long int) igraph_stack_pop(&stack);
      neis = igraph_adjlist_get(adjlist_in_p, actnode);
      nneis = igraph_vector_int_size(neis);
      for (j=0; j<nneis; j++) {
        long int neighbor=(long int) VECTOR(*neis)[j];
	if (nobigint) {
	  tmpscore[neighbor] +=  (tmpscore[actnode]+1)*
	    ((double)(nrgeo[neighbor]))/nrgeo[actnode];
	} else {
	  if (!igraph_biguint_compare_limb(&big_nrgeo[actnode], 0)) {
	    tmpscore[neighbor] = IGRAPH_INFINITY;
	  } else {
	    double div;
	    limb_t shift=1000000000L;
	    IGRAPH_CHECK(igraph_biguint_mul_limb(&T, &big_nrgeo[neighbor], 
						 shift));	  
	    igraph_biguint_div(&D, &R, &T, &big_nrgeo[actnode]);
	    div=igraph_biguint_get(&D) / shift;
	    tmpscore[neighbor] += (tmpscore[actnode]+1) * div;
	  }
	}
      }
      
      if (actnode != source) { VECTOR(*tmpres)[actnode] += tmpscore[actnode]; }

      distance[actnode]=0;
      if (nobigint) { 
	nrgeo[actnode]=0;
      } else {
	igraph_biguint_set_limb(&big_nrgeo[actnode], 0);
      }
      tmpscore[actnode]=0;
      igraph_vector_int_clear(igraph_adjlist_get(adjlist_in_p, actnode));
    }

  } /* for source < no_of_nodes */

  IGRAPH_PROGRESS("Betweenness centrality: ", 100.0, 0);

  /* clean  */
  igraph_Free(distance);
  if (nobigint) {
    igraph_Free(nrgeo); 
  } else {
    igraph_biguint_destroy(&T);
    igraph_biguint_destroy(&R);
    igraph_biguint_destroy(&D);
    IGRAPH_FINALLY_CLEAN(3);
    igraph_i_destroy_biguints(big_nrgeo);
  }
  igraph_Free(tmpscore);
  
  igraph_dqueue_destroy(&q);
  igraph_stack_destroy(&stack);
  IGRAPH_FINALLY_CLEAN(5);

  /* Keep only the requested vertices */
  if (!igraph_vs_is_all(&vids)) { 
    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_VIT_SIZE(vit)));

    for (k=0, IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit);
	 IGRAPH_VIT_NEXT(vit), k++) {
      long int node=IGRAPH_VIT_GET(vit);
      VECTOR(*res)[k] = VECTOR(*tmpres)[node];
    }

    igraph_vit_destroy(&vit);
    igraph_vector_destroy(tmpres);
    IGRAPH_FINALLY_CLEAN(2);
  }     

  /* divide by 2 for undirected graph */
  if (!directed) {
    nneis=igraph_vector_size(res);
    for (j=0; j<nneis; j++) {
      VECTOR(*res)[j] /= 2.0;
    }
  }
  
  igraph_adjlist_destroy(&adjlist_out);
  igraph_adjlist_destroy(&adjlist_in);
  IGRAPH_FINALLY_CLEAN(2);

  return 0;
}

int igraph_i_edge_betweenness_estimate_weighted(const igraph_t *graph, 
					      igraph_vector_t *result,
					      igraph_bool_t directed, 
					      igraph_real_t cutoff,
					      const igraph_vector_t *weights) {

  igraph_real_t minweight;
  igraph_integer_t no_of_nodes=(igraph_integer_t) igraph_vcount(graph);
  igraph_integer_t no_of_edges=(igraph_integer_t) igraph_ecount(graph);
  igraph_2wheap_t Q;
  igraph_inclist_t inclist;
  igraph_inclist_t fathers;
  igraph_neimode_t mode= directed ? IGRAPH_OUT : IGRAPH_ALL;
  igraph_vector_t distance, tmpscore;
  igraph_vector_long_t nrgeo;
  long int source, j;
  int cmp_result;
  const double eps = IGRAPH_SHORTEST_PATH_EPSILON;
  igraph_stack_t S;

  if (igraph_vector_size(weights) != no_of_edges) {
    IGRAPH_ERROR("Weight vector length does not match", IGRAPH_EINVAL);
  }
  minweight = igraph_vector_min(weights);
  if (minweight <= 0) {
    IGRAPH_ERROR("Weight vector must be positive", IGRAPH_EINVAL);
  }
  else if (minweight <= eps) {
    IGRAPH_WARNING("Some weights are smaller than epsilon, calculations may suffer from numerical precision.");
  }
  
  IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, mode));
  IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
  IGRAPH_CHECK(igraph_inclist_init_empty(&fathers, no_of_nodes));
  IGRAPH_FINALLY(igraph_inclist_destroy, &fathers);

  IGRAPH_VECTOR_INIT_FINALLY(&distance, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&tmpscore, no_of_nodes);
  IGRAPH_CHECK(igraph_vector_long_init(&nrgeo, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &nrgeo);

  IGRAPH_CHECK(igraph_2wheap_init(&Q, no_of_nodes));
  IGRAPH_FINALLY(igraph_2wheap_destroy, &Q);
  IGRAPH_CHECK(igraph_stack_init(&S, no_of_nodes));
  IGRAPH_FINALLY(igraph_stack_destroy, &S);

  IGRAPH_CHECK(igraph_vector_resize(result, no_of_edges));
  igraph_vector_null(result);

  for (source=0; source<no_of_nodes; source++) {
    IGRAPH_PROGRESS("Edge betweenness centrality: ", 100.0*source/no_of_nodes, 0);
    IGRAPH_ALLOW_INTERRUPTION();

/*     printf("source: %li\n", source); */
    
    igraph_vector_null(&distance);
    igraph_vector_null(&tmpscore);
    igraph_vector_long_null(&nrgeo);
    
    igraph_2wheap_push_with_index(&Q, source, -1.0);
    VECTOR(distance)[source]=1.0;
    VECTOR(nrgeo)[source]=1;
    
    while (!igraph_2wheap_empty(&Q)) {
      long int minnei=igraph_2wheap_max_index(&Q);
      igraph_real_t mindist=-igraph_2wheap_delete_max(&Q);
      igraph_vector_int_t *neis;
      long int nlen;

      /* printf("SP to %li is final, dist: %g, nrgeo: %li\n", minnei, */
      /* VECTOR(distance)[minnei]-1.0, VECTOR(nrgeo)[minnei]); */
      
      igraph_stack_push(&S, minnei);

      if (cutoff > 0 && VECTOR(distance)[minnei] >= cutoff+1.0) { continue; }

      neis=igraph_inclist_get(&inclist, minnei);
      nlen=igraph_vector_int_size(neis);
      for (j=0; j<nlen; j++) {
	long int edge=(long int) VECTOR(*neis)[j];
	long int to=IGRAPH_OTHER(graph, edge, minnei);
	igraph_real_t altdist=mindist + VECTOR(*weights)[edge];
	igraph_real_t curdist=VECTOR(distance)[to];

	if (curdist == 0) {
	  /* this means curdist is infinity */
	  cmp_result = -1;
	} else {
	  cmp_result = igraph_cmp_epsilon(altdist, curdist, eps);
	}
	
	/* printf("to=%ld, altdist = %lg, curdist = %lg, cmp = %d\n",
	  to, altdist, curdist-1, cmp_result); */
	if (curdist == 0) {
	  /* This is the first finite distance to 'to' */
	  igraph_vector_int_t *v=igraph_inclist_get(&fathers, to);
	  /* printf("Found first path to %li (from %li)\n", to, minnei); */
	  igraph_vector_int_resize(v,1);
	  VECTOR(*v)[0]=edge;
	  VECTOR(nrgeo)[to] = VECTOR(nrgeo)[minnei];
	  VECTOR(distance)[to]=altdist;
	  IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, to, -altdist));
	} else if (cmp_result < 0) {
	  /* This is a shorter path */
	  igraph_vector_int_t *v =igraph_inclist_get(&fathers, to);
	  /* printf("Found a shorter path to %li (from %li)\n", to, minnei); */
	  igraph_vector_int_resize(v,1);
	  VECTOR(*v)[0]=edge;
	  VECTOR(nrgeo)[to] = VECTOR(nrgeo)[minnei];
	  VECTOR(distance)[to] = altdist;
	  IGRAPH_CHECK(igraph_2wheap_modify(&Q, to, -altdist));
	} else if (cmp_result == 0) {
	  igraph_vector_int_t *v=igraph_inclist_get(&fathers, to);
	  /* printf("Found a second SP to %li (from %li)\n", to, minnei); */
	  igraph_vector_int_push_back(v, edge);
	  VECTOR(nrgeo)[to] += VECTOR(nrgeo)[minnei];
	}
      }
	  
    } /* igraph_2wheap_empty(&Q) */

    while (!igraph_stack_empty(&S)) {
      long int w=(long int) igraph_stack_pop(&S);
      igraph_vector_int_t *fatv=igraph_inclist_get(&fathers, w);
      long int fatv_len=igraph_vector_int_size(fatv);
      /* printf("Popping %li.\n", w); */
      for (j=0; j<fatv_len; j++) {
	long int fedge=(long int) VECTOR(*fatv)[j];
	long int neighbor=IGRAPH_OTHER(graph, fedge, w);
	VECTOR(tmpscore)[neighbor] += ((double)VECTOR(nrgeo)[neighbor]) /
	  VECTOR(nrgeo)[w] * (1.0+VECTOR(tmpscore)[w]);
	/* printf("Scoring %li (edge %li)\n", neighbor, fedge); */
	VECTOR(*result)[fedge] += 
	  ((VECTOR(tmpscore)[w]+1) * VECTOR(nrgeo)[neighbor]) / 
	  VECTOR(nrgeo)[w];
      }
      
      VECTOR(tmpscore)[w]=0;
      VECTOR(distance)[w]=0;
      VECTOR(nrgeo)[w]=0;
      igraph_vector_int_clear(fatv);
    }
    
  } /* source < no_of_nodes */

  if (!directed || !igraph_is_directed(graph)) {
    for (j=0; j<no_of_edges; j++) {
      VECTOR(*result)[j] /= 2.0;
    }
  }

  IGRAPH_PROGRESS("Edge betweenness centrality: ", 100.0, 0);

  igraph_stack_destroy(&S);
  igraph_2wheap_destroy(&Q);
  IGRAPH_FINALLY_CLEAN(2);
  
  igraph_inclist_destroy(&inclist);
  igraph_inclist_destroy(&fathers);  
  igraph_vector_destroy(&distance);
  igraph_vector_destroy(&tmpscore);
  igraph_vector_long_destroy(&nrgeo);
  IGRAPH_FINALLY_CLEAN(5);
  
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
 * \param weights An optional weight vector for weighted edge
 *        betweenness. Supply a null pointer here for the unweighted
 *        version. 
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
 * 
 * \example examples/simple/igraph_edge_betweenness.c
 */
int igraph_edge_betweenness(const igraph_t *graph, igraph_vector_t *result,
                            igraph_bool_t directed, 
			    const igraph_vector_t *weights) {
  return igraph_edge_betweenness_estimate(graph, result, directed, -1, 
					  weights);
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
 * \param weights An optional weight vector for weighted
 *        betweenness. Supply a null pointer here for unweighted
 *        betweenness.
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
                                     igraph_bool_t directed, igraph_real_t cutoff,
				     const igraph_vector_t *weights) {
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  long int *distance;
  unsigned long long int *nrgeo;
  double *tmpscore;
  igraph_stack_t stack=IGRAPH_STACK_NULL;
  long int source;
  long int j;

  igraph_inclist_t elist_out, elist_in;
  igraph_inclist_t *elist_out_p, *elist_in_p;
  igraph_vector_int_t *neip;
  long int neino;
  long int i;

  if (weights) { 
    return igraph_i_edge_betweenness_estimate_weighted(graph, result, 
						     directed, cutoff, weights);
  }

  directed=directed && igraph_is_directed(graph);
  if (directed) {
    IGRAPH_CHECK(igraph_inclist_init(graph, &elist_out, IGRAPH_OUT));
    IGRAPH_FINALLY(igraph_inclist_destroy, &elist_out);
    IGRAPH_CHECK(igraph_inclist_init(graph, &elist_in, IGRAPH_IN));
    IGRAPH_FINALLY(igraph_inclist_destroy, &elist_in);
    elist_out_p=&elist_out;
    elist_in_p=&elist_in;
  } else {
    IGRAPH_CHECK(igraph_inclist_init(graph,&elist_out, IGRAPH_ALL));
    IGRAPH_FINALLY(igraph_inclist_destroy, &elist_out);
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

    memset(distance, 0, (size_t) no_of_nodes * sizeof(long int));
    memset(nrgeo, 0, (size_t) no_of_nodes*sizeof(unsigned long long int));
    memset(tmpscore, 0, (size_t) no_of_nodes*sizeof(double));
    igraph_stack_clear(&stack); /* it should be empty anyway... */
    
    IGRAPH_CHECK(igraph_dqueue_push(&q, source));
      
    nrgeo[source]=1;
    distance[source]=0;
    
    while (!igraph_dqueue_empty(&q)) {
      long int actnode=(long int) igraph_dqueue_pop(&q);

      if (cutoff > 0 && distance[actnode] >= cutoff ) continue;

      /* check the neighbors and add to them to the queue if unseen before */
      neip=igraph_inclist_get(elist_out_p, actnode);
      neino=igraph_vector_int_size(neip);
      for (i=0; i<neino; i++) {
	      igraph_integer_t edge=(igraph_integer_t) VECTOR(*neip)[i], from, to;
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
      long int actnode=(long int) igraph_stack_pop(&stack);
      if (distance[actnode]<1) { continue; } /* skip source node */
      
      /* set the temporary score of the friends */
      neip=igraph_inclist_get(elist_in_p, actnode);
      neino=igraph_vector_int_size(neip);
      for (i=0; i<neino; i++) {
	igraph_integer_t from, to;
	long int neighbor;
	igraph_integer_t edgeno=(igraph_integer_t) VECTOR(*neip)[i];
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
    igraph_inclist_destroy(&elist_out);
    igraph_inclist_destroy(&elist_in);
    IGRAPH_FINALLY_CLEAN(2);
  } else {
    igraph_inclist_destroy(&elist_out);
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
 * can be reached from the other vertices). It is defined as
 * the number of vertices minus one divided by the sum of the
 * lengths of all geodesics from/to the given vertex.
 *
 * </para><para>
 * If the graph is not connected, and there is no path between two
 * vertices, the number of vertices is used instead the length of the
 * geodesic. This is longer than the longest possible geodesic in case
 * of unweighted graphs, but may not be so in weighted graphs, so it is
 * best not to use this function on weighted graphs.
 * 
 * </para><para>
 * If the graph has a single vertex only, the closeness centrality of
 * that single vertex will be NaN (because we are essentially dividing
 * zero with zero).
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
 * \param weights An optional vector containing edge weights for
 *        weighted closeness. Supply a null pointer here for
 *        traditional, unweighted closeness.
 * \param normalized Boolean, whether to normalize results by multiplying
 *        by the number of vertices minus one.
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
                     const igraph_vs_t vids, igraph_neimode_t mode, 
		     const igraph_vector_t *weights,
		     igraph_bool_t normalized) {
  return igraph_closeness_estimate(graph, res, vids, mode, -1, weights,
				   normalized);
}

int igraph_i_closeness_estimate_weighted(const igraph_t *graph, 
				       igraph_vector_t *res, 
				       const igraph_vs_t vids, 
				       igraph_neimode_t mode,
				       igraph_real_t cutoff,
				       const igraph_vector_t *weights,
				       igraph_bool_t normalized) {

  /* See igraph_shortest_paths_dijkstra() for the implementation 
     details and the dirty tricks. */

  igraph_real_t minweight;
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  
  igraph_2wheap_t Q;
  igraph_vit_t vit;
  long int nodes_to_calc;
  
  igraph_lazy_inclist_t inclist;
  long int i, j;
  
  igraph_vector_t dist;
  igraph_vector_long_t which;
  long int nodes_reached;

  int cmp_result;
  const double eps = IGRAPH_SHORTEST_PATH_EPSILON;
  igraph_real_t mindist;

  igraph_bool_t warning_shown = 0;
  
  if (igraph_vector_size(weights) != no_of_edges) {
    IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
  }
  
  minweight = igraph_vector_min(weights);
  if (minweight <= 0) {
    IGRAPH_ERROR("Weight vector must be positive", IGRAPH_EINVAL);
  }
  else if (minweight <= eps) {
    IGRAPH_WARNING("Some weights are smaller than epsilon, calculations may suffer from numerical precision.");
  }
  
  IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);
  
  nodes_to_calc=IGRAPH_VIT_SIZE(vit);
  
  IGRAPH_CHECK(igraph_2wheap_init(&Q, no_of_nodes));
  IGRAPH_FINALLY(igraph_2wheap_destroy, &Q);
  IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, mode));
  IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

  IGRAPH_VECTOR_INIT_FINALLY(&dist, no_of_nodes);
  IGRAPH_CHECK(igraph_vector_long_init(&which, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &which);

  IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));
  igraph_vector_null(res);

  for (i=0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
    
    long int source=IGRAPH_VIT_GET(vit);
    igraph_2wheap_clear(&Q);
    igraph_2wheap_push_with_index(&Q, source, -1.0);
    VECTOR(which)[source] = i+1;
    VECTOR(dist)[source] = 1.0;     /* actual distance is zero but we need to store distance + 1 */
    nodes_reached=0;

    while (!igraph_2wheap_empty(&Q)) {
      igraph_integer_t minnei=(igraph_integer_t) igraph_2wheap_max_index(&Q);
      /* Now check all neighbors of minnei for a shorter path */
      igraph_vector_t *neis=igraph_lazy_inclist_get(&inclist, minnei);
      long int nlen=igraph_vector_size(neis);

      mindist=-igraph_2wheap_delete_max(&Q);

      VECTOR(*res)[i] += (mindist - 1.0);
      nodes_reached++;

      if (cutoff>0 && mindist>=cutoff+1.0) continue;    /* NOT break!!! */

      for (j=0; j<nlen; j++) {
	      long int edge=(long int) VECTOR(*neis)[j];
	      long int to=IGRAPH_OTHER(graph, edge, minnei);
	      igraph_real_t altdist=mindist+VECTOR(*weights)[edge];
	      igraph_real_t curdist=VECTOR(dist)[to];
	      if (curdist == 0) {
	        /* this means curdist is infinity */
	        cmp_result = -1;
	      } else {
	        cmp_result = igraph_cmp_epsilon(altdist, curdist, eps);
	      }

	      if (VECTOR(which)[to] != i+1) {
	        /* First non-infinite distance */
	        VECTOR(which)[to]=i+1;
	        VECTOR(dist)[to]=altdist;
	        IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, to, -altdist));
	      } else if (cmp_result < 0) {
	        /* This is a shorter path */
	        VECTOR(dist)[to]=altdist;
	        IGRAPH_CHECK(igraph_2wheap_modify(&Q, to, -altdist));
	      }
      }

    } /* !igraph_2wheap_empty(&Q) */

    /* using igraph_real_t here instead of igraph_integer_t to avoid overflow */
    VECTOR(*res)[i] += ((igraph_real_t)no_of_nodes * (no_of_nodes-nodes_reached));
    VECTOR(*res)[i] = (no_of_nodes-1) / VECTOR(*res)[i];

    if (((cutoff > 0 && mindist < cutoff+1.0) || (cutoff <= 0)) &&
        nodes_reached < no_of_nodes && !warning_shown) {
      IGRAPH_WARNING("closeness centrality is not well-defined for disconnected graphs");
      warning_shown = 1;
    }
  } /* !IGRAPH_VIT_END(vit) */

  if (!normalized) {
    for (i=0; i<nodes_to_calc; i++) {
      VECTOR(*res)[i] /= (no_of_nodes-1);
    }
  }

  igraph_vector_long_destroy(&which);
  igraph_vector_destroy(&dist);
  igraph_lazy_inclist_destroy(&inclist);
  igraph_2wheap_destroy(&Q);
  igraph_vit_destroy(&vit);
  IGRAPH_FINALLY_CLEAN(5);

  return 0;
}

/**
 * \ingroup structural
 * \function igraph_closeness_estimate
 * \brief Closeness centrality estimations for some vertices.
 *
 * </para><para>
 * The closeness centrality of a vertex measures how easily other
 * vertices can be reached from it (or the other way: how easily it
 * can be reached from the other vertices). It is defined as
 * the number of vertices minus one divided by the sum of the
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
 * \param weights An optional vector containing edge weights for
 *        weighted closeness. Supply a null pointer here for
 *        traditional, unweighted closeness.
 * \param normalized Boolean, whether to normalize results by multiplying
 *        by the number of vertices minus one.
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
                              igraph_real_t cutoff,
			      const igraph_vector_t *weights,
			      igraph_bool_t normalized) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_vector_t already_counted;
  igraph_vector_int_t *neis;
  long int i, j;
  long int nodes_reached;
  long int actdist;
  igraph_adjlist_t allneis;

  igraph_dqueue_t q;
  
  long int nodes_to_calc;
  igraph_vit_t vit;

  igraph_bool_t warning_shown = 0;
  
  if (weights) { 
    return igraph_i_closeness_estimate_weighted(graph, res, vids, mode, cutoff,
						weights, normalized);
  }

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
    igraph_dqueue_clear(&q);
    IGRAPH_CHECK(igraph_dqueue_push(&q, IGRAPH_VIT_GET(vit)));
    IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
    nodes_reached=1;
    VECTOR(already_counted)[(long int)IGRAPH_VIT_GET(vit)]=i+1;

    IGRAPH_PROGRESS("Closeness: ", 100.0*i/no_of_nodes, NULL);
    IGRAPH_ALLOW_INTERRUPTION();
    
    while (!igraph_dqueue_empty(&q)) {
      long int act=(long int) igraph_dqueue_pop(&q);
      actdist=(long int) igraph_dqueue_pop(&q);
      
      VECTOR(*res)[i] += actdist;

      if (cutoff > 0 && actdist >= cutoff) continue;   /* NOT break!!! */

      /* check the neighbors */
      neis=igraph_adjlist_get(&allneis, act);
      for (j=0; j<igraph_vector_int_size(neis); j++) {
        long int neighbor=(long int) VECTOR(*neis)[j];
        if (VECTOR(already_counted)[neighbor] == i+1) { continue; }
        VECTOR(already_counted)[neighbor] = i+1;
        nodes_reached++;
        IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
        IGRAPH_CHECK(igraph_dqueue_push(&q, actdist+1));
      }
    }

    /* using igraph_real_t here instead of igraph_integer_t to avoid overflow */
    VECTOR(*res)[i] += ((igraph_real_t)no_of_nodes * (no_of_nodes-nodes_reached));
    VECTOR(*res)[i] = (no_of_nodes-1) / VECTOR(*res)[i];

    if (((cutoff > 0 && actdist < cutoff) || cutoff <= 0) &&
        no_of_nodes > nodes_reached && !warning_shown) {
      IGRAPH_WARNING("closeness centrality is not well-defined for disconnected graphs");
      warning_shown = 1;
    }
  }

  if (!normalized) {
    for (i=0; i<nodes_to_calc; i++) {
      VECTOR(*res)[i] /= (no_of_nodes-1);
    }
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

/**
 * \function igraph_centralization
 * Calculate the centralization score from the node level scores
 * 
 * For a centrality score defined on the vertices of a graph, it is
 * possible to define a graph level centralization index, by
 * calculating the sum of the deviation from the maximum centrality
 * score. Consequently, the higher the centralization index of the
 * graph, the more centralized the structure is.
 * 
 * </para><para>In order to make graphs of different sizes comparable,
 * the centralization index is usually normalized to a number between
 * zero and one, by dividing the (unnormalized) centralization score
 * of the most centralized structure with the same number of vertices.
 * 
 * </para><para>For most centrality indices the most centralized
 * structure is the star graph, a single center connected to all other
 * nodes in the network. There are some variation depending on whether
 * the graph is directed or not, whether loop edges are allowed, etc. 
 * 
 * </para><para>
 * This function simply calculates the graph level index, if the node
 * level scores and the theoretical maximum are given. It is called by
 * all the measure-specific centralization functions.
 * 
 * \param scores A vector containing the node-level centrality
 *     scores.
 * \param theoretical_max The graph level centrality score of the most
 *     centralized graph with the same number of vertices. Only used
 *     if \c normalized set to true.
 * \param normalized Boolean, whether to normalize the centralization
 *     by dividing the supplied theoretical maximum.
 * \return The graph level index.
 * 
 * \sa \ref igraph_centralization_degree(), \ref
 * igraph_centralization_betweenness(), \ref
 * igraph_centralization_closeness(), and \ref
 * igraph_centralization_eigenvector_centrality() for specific
 * centralization functions.
 * 
 * Time complexity: O(n), the length of the score vector.
 * 
 * \example examples/simple/centralization.c
 */

igraph_real_t igraph_centralization(const igraph_vector_t *scores,
				    igraph_real_t theoretical_max,
				    igraph_bool_t normalized) {
  
  long int no_of_nodes=igraph_vector_size(scores);
  igraph_real_t maxscore=0.0;
  igraph_real_t cent=0.0;
  
  if (no_of_nodes != 0) {
    maxscore = igraph_vector_max(scores);
    cent = no_of_nodes * maxscore - igraph_vector_sum(scores);
    if (normalized) { cent = cent/theoretical_max; }
  } else {
    cent = IGRAPH_NAN;
  }

  return cent;
}

/**
 * \function igraph_centralization_degree
 * Calculate vertex degree and graph centralization
 * 
 * This function calculates the degree of the vertices by passing its
 * arguments to \ref igraph_degree(); and it calculates the graph
 * level centralization index based on the results by calling \ref
 * igraph_centralization().
 * \param graph The input graph.
 * \param res A vector if you need the node-level degree scores, or a
 *     null pointer otherwise.
 * \param mode Constant the specifies the type of degree for directed
 *     graphs. Possible values: \c IGRAPH_IN, \c IGRAPH_OUT and \c
 *     IGRAPH_ALL. This argument is ignored for undirected graphs.
 * \param loops Boolean, whether to consider loop edges when
 *     calculating the degree (and the centralization).
 * \param centralization Pointer to a real number, the centralization
 *     score is placed here.
 * \param theoretical_max Pointer to real number or a null pointer. If
 *     not a null pointer, then the theoretical maximum graph
 *     centrality score for a graph with the same number vertices is
 *     stored here.
 * \param normalized Boolean, whether to calculate a normalized
 *     centralization score. See \ref igraph_centralization() for how
 *     the normalization is done.
 * \return Error code.
 * 
 * \sa \ref igraph_centralization(), \ref igraph_degree().
 * 
 * Time complexity: the complexity of \ref igraph_degree() plus O(n),
 * the number of vertices queried, for calculating the centralization
 * score.
 */

int igraph_centralization_degree(const igraph_t *graph, igraph_vector_t *res, 
				 igraph_neimode_t mode, igraph_bool_t loops,
				 igraph_real_t *centralization,
				 igraph_real_t *theoretical_max,
				 igraph_bool_t normalized) {
  
  igraph_vector_t myscores;
  igraph_vector_t *scores=res;
  igraph_real_t *tmax=theoretical_max, mytmax;

  if (!tmax) { tmax=&mytmax; }

  if (!res) {
    scores=&myscores;
    IGRAPH_VECTOR_INIT_FINALLY(scores, 0);
  }
  
  IGRAPH_CHECK(igraph_degree(graph, scores, igraph_vss_all(), mode, loops));
  
  IGRAPH_CHECK(igraph_centralization_degree_tmax(graph, 0, mode, loops, 
						 tmax));

  *centralization = igraph_centralization(scores, *tmax, normalized);
  
  if (!res) {
    igraph_vector_destroy(scores);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  return 0;
}

/** 
 * \function igraph_centralization_degree_tmax
 * Theoretical maximum for graph centralization based on degree
 * 
 * This function returns the theoretical maximum graph centrality
 * based on vertex degree. 
 * 
 * </para><para>
 * There are two ways to call this function, the first is to supply a
 * graph as the <code>graph</code> argument, and then the number of
 * vertices is taken from this object, and its directedness is
 * considered as well. The <code>nodes</code> argument is ignored in
 * this case. The <code>mode</code> argument is also ignored if the
 * supplied graph is undirected.
 * 
 * </para><para>
 * The other way is to supply a null pointer as the <code>graph</code>
 * argument. In this case the <code>nodes</code> and <code>mode</code>
 * arguments are considered.
 * 
 * </para><para>
 * The most centralized structure is the star. More specifically, for
 * undirected graphs it is the star, for directed graphs it is the
 * in-star or the out-star.
 * \param graph A graph object or a null pointer, see the description
 *     above.
 * \param nodes The number of nodes. This is ignored if the
 *     <code>graph</code> argument is not a null pointer.
 * \param mode Constant, whether the calculation is based on in-degree
 *     (<code>IGRAPH_IN</code>), out-degree (<code>IGRAPH_OUT</code>)
 *     or total degree (<code>IGRAPH_ALL</code>). This is ignored if
 *     the <code>graph</code> argument is not a null pointer and the
 *     given graph is undirected. 
 * \param loops Boolean scalar, whether to consider loop edges in the
 *     calculation. 
 * \param res Pointer to a real variable, the result is stored here.
 * \return Error code.
 * 
 * Time complexity: O(1).
 * 
 * \sa \ref igraph_centralization_degree() and \ref
 * igraph_centralization().
 */

int igraph_centralization_degree_tmax(const igraph_t *graph, 
				      igraph_integer_t nodes,
				      igraph_neimode_t mode,
				      igraph_bool_t loops,
				      igraph_real_t *res) {

  igraph_bool_t directed=mode != IGRAPH_ALL;
  igraph_real_t real_nodes;

  if (graph) {
    directed=igraph_is_directed(graph);
    nodes=igraph_vcount(graph);
  }

  real_nodes = nodes;    /* implicit cast to igraph_real_t */

  if (directed) {
    switch (mode) {
    case IGRAPH_IN:
    case IGRAPH_OUT:
      if (!loops) {
	*res = (real_nodes-1) * (real_nodes-1);
      } else {
	*res = (real_nodes-1) * real_nodes;
      }
      break;
    case IGRAPH_ALL:
      if (!loops) {
	*res = 2 * (real_nodes-1) * (real_nodes-2);
      } else {
	*res = 2 * (real_nodes-1) * (real_nodes-1);
      }
      break;
    }
  } else {
    if (!loops) {
      *res = (real_nodes-1) * (real_nodes-2);
    } else {
      *res = (real_nodes-1) * real_nodes;
    }
  }

  return 0;
}

/**
 * \function igraph_centralization_betweenness
 * Calculate vertex betweenness and graph centralization
 * 
 * This function calculates the betweenness centrality of the vertices
 * by passing its arguments to \ref igraph_betweenness(); and it
 * calculates the graph level centralization index based on the
 * results by calling \ref igraph_centralization().
 * \param graph The input graph.
 * \param res A vector if you need the node-level betweenness scores, or a
 *     null pointer otherwise.
 * \param directed Boolean, whether to consider directed paths when
 *     calculating betweenness.
 * \param nobigint Logical, if true, then we don't use big integers
 *        for the calculation, setting this to zero (=false) should
 *        work for most graphs. It is currently ignored for weighted
 *        graphs.
 * \param centralization Pointer to a real number, the centralization
 *     score is placed here.
 * \param theoretical_max Pointer to real number or a null pointer. If
 *     not a null pointer, then the theoretical maximum graph
 *     centrality score for a graph with the same number vertices is
 *     stored here.
 * \param normalized Boolean, whether to calculate a normalized
 *     centralization score. See \ref igraph_centralization() for how
 *     the normalization is done.
 * \return Error code.
 * 
 * \sa \ref igraph_centralization(), \ref igraph_betweenness().
 * 
 * Time complexity: the complexity of \ref igraph_betweenness() plus
 * O(n), the number of vertices queried, for calculating the
 * centralization score.
 */

int igraph_centralization_betweenness(const igraph_t *graph, 
				      igraph_vector_t *res,
				      igraph_bool_t directed,
				      igraph_bool_t nobigint,
				      igraph_real_t *centralization,
				      igraph_real_t *theoretical_max,
				      igraph_bool_t normalized) {
  
  igraph_vector_t myscores;
  igraph_vector_t *scores=res;
  igraph_real_t *tmax=theoretical_max, mytmax;

  if (!tmax) { tmax=&mytmax; }

  if (!res) {
    scores=&myscores;
    IGRAPH_VECTOR_INIT_FINALLY(scores, 0);
  }
  
  IGRAPH_CHECK(igraph_betweenness(graph, scores, igraph_vss_all(), directed, 
				  /*weights=*/ 0, nobigint));
  
  IGRAPH_CHECK(igraph_centralization_betweenness_tmax(graph, 0, directed, 
						      tmax));
  
  *centralization = igraph_centralization(scores, *tmax, normalized);
  
  if (!res) {
    igraph_vector_destroy(scores);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  return 0;
}

/** 
 * \function igraph_centralization_betweenness_tmax
 * Theoretical maximum for graph centralization based on betweenness
 * 
 * This function returns the theoretical maximum graph centrality
 * based on vertex betweenness. 
 * 
 * </para><para>
 * There are two ways to call this function, the first is to supply a
 * graph as the <code>graph</code> argument, and then the number of
 * vertices is taken from this object, and its directedness is
 * considered as well. The <code>nodes</code> argument is ignored in
 * this case. The <code>directed</code> argument is also ignored if the
 * supplied graph is undirected.
 * 
 * </para><para>
 * The other way is to supply a null pointer as the <code>graph</code>
 * argument. In this case the <code>nodes</code> and <code>directed</code>
 * arguments are considered.
 * 
 * </para><para>
 * The most centralized structure is the star.
 * \param graph A graph object or a null pointer, see the description
 *     above.
 * \param nodes The number of nodes. This is ignored if the
 *     <code>graph</code> argument is not a null pointer.
 * \param directed Boolean scalar, whether to use directed paths in
 *     the betweenness calculation. This argument is ignored if
 *     <code>graph</code> is not a null pointer and it is undirected.
 * \param res Pointer to a real variable, the result is stored here.
 * \return Error code.
 * 
 * Time complexity: O(1).
 * 
 * \sa \ref igraph_centralization_betweenness() and \ref
 * igraph_centralization().
 */

int igraph_centralization_betweenness_tmax(const igraph_t *graph, 
					   igraph_integer_t nodes,
					   igraph_bool_t directed,
					   igraph_real_t *res) {
  igraph_real_t real_nodes;

  if (graph) { 
    directed=directed && igraph_is_directed(graph); 
    nodes=igraph_vcount(graph);
  }

  real_nodes = nodes;    /* implicit cast to igraph_real_t */

  if (directed) {
    *res = (real_nodes-1) * (real_nodes-1) * (real_nodes-2);
  } else {
    *res = (real_nodes-1) * (real_nodes-1) * (real_nodes-2) / 2.0;
  }
  
  return 0;
}

/**
 * \function igraph_centralization_closeness
 * Calculate vertex closeness and graph centralization
 * 
 * This function calculates the closeness centrality of the vertices
 * by passing its arguments to \ref igraph_closeness(); and it
 * calculates the graph level centralization index based on the
 * results by calling \ref igraph_centralization().
 * \param graph The input graph.
 * \param res A vector if you need the node-level closeness scores, or a
 *     null pointer otherwise.
 * \param mode Constant the specifies the type of closeness for directed
 *     graphs. Possible values: \c IGRAPH_IN, \c IGRAPH_OUT and \c
 *     IGRAPH_ALL. This argument is ignored for undirected graphs. See
 *     \ref igraph_closeness() argument with the same name for more.
 * \param centralization Pointer to a real number, the centralization
 *     score is placed here.
 * \param theoretical_max Pointer to real number or a null pointer. If
 *     not a null pointer, then the theoretical maximum graph
 *     centrality score for a graph with the same number vertices is
 *     stored here.
 * \param normalized Boolean, whether to calculate a normalized
 *     centralization score. See \ref igraph_centralization() for how
 *     the normalization is done.
 * \return Error code.
 * 
 * \sa \ref igraph_centralization(), \ref igraph_closeness().
 * 
 * Time complexity: the complexity of \ref igraph_closeness() plus
 * O(n), the number of vertices queried, for calculating the
 * centralization score.
 */

int igraph_centralization_closeness(const igraph_t *graph, 
				    igraph_vector_t *res, 
				    igraph_neimode_t mode, 
				    igraph_real_t *centralization,
				    igraph_real_t *theoretical_max,
				    igraph_bool_t normalized) {

  igraph_vector_t myscores;
  igraph_vector_t *scores=res;
  igraph_real_t *tmax=theoretical_max, mytmax;

  if (!tmax) { tmax=&mytmax; }    

  if (!res) {
    scores=&myscores;
    IGRAPH_VECTOR_INIT_FINALLY(scores, 0);
  }
  
  IGRAPH_CHECK(igraph_closeness(graph, scores, igraph_vss_all(), mode, 
				/*weights=*/ 0, /*normalize=*/ 1));

  IGRAPH_CHECK(igraph_centralization_closeness_tmax(graph, 0, mode, 
						    tmax));

  *centralization = igraph_centralization(scores, *tmax, normalized);
  
  if (!res) {
    igraph_vector_destroy(scores);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  return 0;
}

/** 
 * \function igraph_centralization_closeness_tmax
 * Theoretical maximum for graph centralization based on closeness
 * 
 * This function returns the theoretical maximum graph centrality
 * based on vertex closeness. 
 * 
 * </para><para>
 * There are two ways to call this function, the first is to supply a
 * graph as the <code>graph</code> argument, and then the number of
 * vertices is taken from this object, and its directedness is
 * considered as well. The <code>nodes</code> argument is ignored in
 * this case. The <code>mode</code> argument is also ignored if the
 * supplied graph is undirected.
 * 
 * </para><para>
 * The other way is to supply a null pointer as the <code>graph</code>
 * argument. In this case the <code>nodes</code> and <code>mode</code>
 * arguments are considered.
 * 
 * </para><para>
 * The most centralized structure is the star.
 * \param graph A graph object or a null pointer, see the description
 *     above.
 * \param nodes The number of nodes. This is ignored if the
 *     <code>graph</code> argument is not a null pointer.
 * \param mode Constant, specifies what kinf of distances to consider
 *     to calculate closeness. See the <code>mode</code> argument of
 *     \ref igraph_closeness() for details. This argument is ignored
 *     if <code>graph</code> is not a null pointer and it is
 *     undirected.
 * \param res Pointer to a real variable, the result is stored here.
 * \return Error code.
 * 
 * Time complexity: O(1).
 * 
 * \sa \ref igraph_centralization_closeness() and \ref
 * igraph_centralization().
 */

int igraph_centralization_closeness_tmax(const igraph_t *graph,
					 igraph_integer_t nodes,
					 igraph_neimode_t mode,
					 igraph_real_t *res) {
  igraph_real_t real_nodes;

  if (graph) {
    nodes=igraph_vcount(graph); 
    if (!igraph_is_directed(graph)) { mode=IGRAPH_ALL; }
  }

  real_nodes = nodes;    /* implicit cast to igraph_real_t */

  if (mode != IGRAPH_ALL) {
    *res = (real_nodes-1) * (1.0-1.0/real_nodes);
  } else {
    *res = (real_nodes-1) * (real_nodes-2) / (2.0*real_nodes-3);
  }
  
  return 0;
}

/**
 * \function igraph_centralization_eigenvector_centrality
 * Calculate eigenvector centrality scores and graph centralization
 * 
 * This function calculates the eigenvector centrality of the vertices
 * by passing its arguments to \ref igraph_eigenvector_centrality);
 * and it calculates the graph level centralization index based on the
 * results by calling \ref igraph_centralization().
 * \param graph The input graph.
 * \param vector A vector if you need the node-level eigenvector
 *      centrality scores, or a null pointer otherwise.
 * \param value If not a null pointer, then the leading eigenvalue is
 *      stored here.
 * \param scale If not zero then the result will be scaled, such that
 *     the absolute value of the maximum centrality is one.
 * \param options Options to ARPACK. See \ref igraph_arpack_options_t
 *    for details. Note that the function overwrites the
 *    <code>n</code> (number of vertices) parameter and 
 *    it always starts the calculation from a non-random vector
 *    calculated based on the degree of the vertices.
 * \param centralization Pointer to a real number, the centralization
 *     score is placed here.
 * \param theoretical_max Pointer to real number or a null pointer. If
 *     not a null pointer, then the theoretical maximum graph
 *     centrality score for a graph with the same number vertices is
 *     stored here.
 * \param normalized Boolean, whether to calculate a normalized
 *     centralization score. See \ref igraph_centralization() for how
 *     the normalization is done.
 * \return Error code.
 * 
 * \sa \ref igraph_centralization(), \ref igraph_eigenvector_centrality().
 * 
 * Time complexity: the complexity of \ref
 * igraph_eigenvector_centrality() plus O(|V|), the number of vertices
 * for the calculating the centralization.
 */

int igraph_centralization_eigenvector_centrality(
					 const igraph_t *graph,
					 igraph_vector_t *vector,
					 igraph_real_t *value,
					 igraph_bool_t directed,
					 igraph_bool_t scale,
					 igraph_arpack_options_t *options,
					 igraph_real_t *centralization,
					 igraph_real_t *theoretical_max,
					 igraph_bool_t normalized) {
  
  igraph_vector_t myscores;
  igraph_vector_t *scores=vector;
  igraph_real_t realvalue, *myvalue=value;
  igraph_real_t *tmax=theoretical_max, mytmax;

  if (!tmax) { tmax=&mytmax; }

  if (!vector) {
    scores=&myscores;
    IGRAPH_VECTOR_INIT_FINALLY(scores, 0);
  }
  if (!value) {
    myvalue=&realvalue;
  }
  
  IGRAPH_CHECK(igraph_eigenvector_centrality(graph, scores, myvalue, directed,
					     scale, /*weights=*/ 0, 
					     options));

  IGRAPH_CHECK(igraph_centralization_eigenvector_centrality_tmax(
						 graph, 0, directed, 
						 scale, 
						 tmax));

  *centralization = igraph_centralization(scores, *tmax, normalized);
  
  if (!vector) {
    igraph_vector_destroy(scores);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  return 0;
}

/**
 * \function igraph_centralization_eigenvector_centrality_tmax
 * Theoretical maximum centralization for eigenvector centrality
 * 
 * This function returns the theoretical maximum graph centrality
 * based on vertex eigenvector centrality. 
 * 
 * </para><para>
 * There are two ways to call this function, the first is to supply a
 * graph as the <code>graph</code> argument, and then the number of
 * vertices is taken from this object, and its directedness is
 * considered as well. The <code>nodes</code> argument is ignored in
 * this case. The <code>directed</code> argument is also ignored if the
 * supplied graph is undirected.
 * 
 * </para><para>
 * The other way is to supply a null pointer as the <code>graph</code>
 * argument. In this case the <code>nodes</code> and <code>directed</code>
 * arguments are considered.
 * 
 * </para><para>
 * The most centralized directed structure is the in-star. The most
 * centralized undirected structure is the graph with a single edge.
 * \param graph A graph object or a null pointer, see the description
 *     above.
 * \param nodes The number of nodes. This is ignored if the
 *     <code>graph</code> argument is not a null pointer.
 * \param directed Boolean scalar, whether to consider edge
 *     directions. This argument is ignored if
 *     <code>graph</code> is not a null pointer and it is undirected.
 * \param scale Whether to rescale the node-level centrality scores to
 *     have a maximum of one.
 * \param res Pointer to a real variable, the result is stored here.
 * \return Error code.
 * 
 * Time complexity: O(1).
 * 
 * \sa \ref igraph_centralization_closeness() and \ref
 * igraph_centralization().
 */

int igraph_centralization_eigenvector_centrality_tmax(
					 const igraph_t *graph,
					 igraph_integer_t nodes,
					 igraph_bool_t directed,
					 igraph_bool_t scale, 
					 igraph_real_t *res) {

  if (graph) {
    nodes=igraph_vcount(graph);
    directed=directed && igraph_is_directed(graph);
  }
  
  if (directed) {
    *res = nodes - 1; 
  } else {
    if (scale) { 
      *res = nodes - 2;
    } else {
      *res = (nodes-2.0) / M_SQRT2;
    }
  }

  return 0;
}
