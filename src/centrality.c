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

int igraph_eigenvector_centrality(const igraph_t *graph, igraph_vector_t *vector,
				  igraph_real_t *value, igraph_bool_t scale,
				  const igraph_vector_t *weights,
				  igraph_arpack_options_t *options) {
  
  igraph_vector_t values;
  igraph_matrix_t vectors;
  
  options->n=igraph_vcount(graph);

  if (weights && igraph_vector_size(weights) != igraph_ecount(graph)) {
    IGRAPH_ERROR("Invalid length of weights vector when calculating "
		 "eigenvector centrality", IGRAPH_EINVAL);
  }

  if (weights && igraph_is_directed(graph)) {
    IGRAPH_WARNING("Weighted directed graph in eigenvector centrality");
  }

  IGRAPH_VECTOR_INIT_FINALLY(&values, 0);
  IGRAPH_MATRIX_INIT_FINALLY(&vectors, 0, 0);

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
  
  options->n=igraph_vcount(graph);
  
  IGRAPH_VECTOR_INIT_FINALLY(&values, 0);
  IGRAPH_MATRIX_INIT_FINALLY(&vectors, 0, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&tmp, options->n);
  
  if (inout==0) {
    inadjlist=&myinadjlist; 
    outadjlist=&myoutadjlist;
  } else if (inout==1) {
    inadjlist=&myoutadjlist;
    outadjlist=&myinadjlist;
  } else {
    /* This should not happen */
    IGRAPH_ERROR("Invalid 'inout' argument, plese do not call "
		 "this funtion directly", IGRAPH_FAILURE);
  }

  IGRAPH_CHECK(igraph_adjlist_init(graph, &myinadjlist, IGRAPH_IN));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &myinadjlist);
  IGRAPH_CHECK(igraph_adjlist_init(graph, &myoutadjlist, IGRAPH_OUT));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &myoutadjlist);

  extra.in=inadjlist; extra.out=outadjlist; extra.tmp=&tmp;

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

int igraph_hub_score(const igraph_t *graph, igraph_vector_t *vector,
		     igraph_real_t *value, igraph_bool_t scale,
		     igraph_arpack_options_t *options) {

  return igraph_i_kleinberg(graph, vector, value, scale, options, 0);
}
			    
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
  
  options->n = igraph_vcount(graph);
  
  if (directed && igraph_is_directed(graph)) {
    IGRAPH_ERROR("Directed graphs are not supported yet, use 'directed=0'",
		 IGRAPH_UNIMPLEMENTED);
  }
  if (weights && igraph_vector_size(weights) != igraph_ecount(graph))
  {
    IGRAPH_ERROR("Invalid length of weights vector when calculating "
		 "eigenvector centrality", IGRAPH_EINVAL);
  }
  
  IGRAPH_MATRIX_INIT_FINALLY(&values, 0, 0);
  IGRAPH_MATRIX_INIT_FINALLY(&vectors, 0, 0);

  if (directed) { dirmode=IGRAPH_IN; } else { dirmode=IGRAPH_ALL; }

  IGRAPH_VECTOR_INIT_FINALLY(&outdegree, options->n);
  IGRAPH_VECTOR_INIT_FINALLY(&tmp, options->n);

  if (!weights) {
    
    igraph_adjlist_t adjlist;
    igraph_i_pagerank_data_t data = { graph, &adjlist, damping,
				      &outdegree, &tmp };

    IGRAPH_CHECK(igraph_degree(graph, &outdegree, igraph_vss_all(),
			       IGRAPH_OUT, /*loops=*/ 0));
    /* Avoid division by zero */
    for (i=0; i<options->n; i++) {
      if (VECTOR(outdegree)[i]==0) {
	VECTOR(outdegree)[i]=1;
      }
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
      VECTOR(outdegree)[to]   += weight;
    }
    
    IGRAPH_CHECK(igraph_arpack_rnsolve(igraph_i_pagerank2,
				       &data, options, 0, &values, &vectors));
    
    igraph_adjedgelist_destroy(&adjedgelist);
    IGRAPH_FINALLY_CLEAN(1);
  }

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
  long int *nrgeo;
  double *tmpscore;
  igraph_stack_t stack=IGRAPH_STACK_NULL;
  long int source;
  long int j, k;
  igraph_vector_t tmp=IGRAPH_VECTOR_NULL;
  igraph_integer_t modein, modeout;
  igraph_vit_t vit;

  IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);

  if (directed) 
    { modeout=IGRAPH_OUT; modein=IGRAPH_IN; } 
  else 
    { modeout=modein=IGRAPH_ALL; }

  distance=igraph_Calloc(no_of_nodes, long int);
  if (distance==0) {
    IGRAPH_ERROR("betweenness failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, distance);
  nrgeo=igraph_Calloc(no_of_nodes, long int);
  if (nrgeo==0) {
    IGRAPH_ERROR("betweenness failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, nrgeo);
  tmpscore=igraph_Calloc(no_of_nodes, double);
  if (tmpscore==0) {
    IGRAPH_ERROR("betweenness failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, tmpscore);

  IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
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
    memset(nrgeo, 0, no_of_nodes*sizeof(long int));
    memset(tmpscore, 0, no_of_nodes*sizeof(double));
    igraph_stack_clear(&stack); /* it should be empty anyway... */
    
    IGRAPH_CHECK(igraph_dqueue_push(&q, source));
    nrgeo[source]=1;
    distance[source]=0;
    
    while (!igraph_dqueue_empty(&q)) {
      long int actnode=igraph_dqueue_pop(&q);

      if (cutoff > 0 && distance[actnode] >= cutoff) continue;
       
      IGRAPH_CHECK(igraph_neighbors(graph, &tmp, actnode, modeout));
      for (j=0; j<igraph_vector_size(&tmp); j++) {
	long int neighbor=VECTOR(tmp)[j];
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
      IGRAPH_CHECK(igraph_neighbors(graph, &tmp, actnode, modein));
      for (j=0; j<igraph_vector_size(&tmp); j++) {
	long int neighbor=VECTOR(tmp)[j];
	if (distance[neighbor]==distance[actnode]-1 &&
	    nrgeo[neighbor] != 0) {
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
  if (!directed || !igraph_is_directed(graph)) {
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
  igraph_vector_destroy(&tmp);
  igraph_vit_destroy(&vit);
  IGRAPH_FINALLY_CLEAN(7);

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
  long int *nrgeo;
  double *tmpscore;
  igraph_stack_t stack=IGRAPH_STACK_NULL;
  long int source;
  long int j;

  igraph_adjedgelist_t elist_out, elist_in;
  igraph_adjedgelist_t *elist_out_p, *elist_in_p;
  igraph_vector_t *neip;
  long int neino;
  long int i;
  igraph_integer_t modein, modeout;

  directed=directed && igraph_is_directed(graph);
  if (directed) {
    modeout=IGRAPH_OUT;
    modein=IGRAPH_IN;
    IGRAPH_CHECK(igraph_adjedgelist_init(graph, &elist_out, IGRAPH_OUT));
    IGRAPH_FINALLY(igraph_adjedgelist_destroy, &elist_out);
    IGRAPH_CHECK(igraph_adjedgelist_init(graph, &elist_in, IGRAPH_IN));
    IGRAPH_FINALLY(igraph_adjedgelist_destroy, &elist_in);
    elist_out_p=&elist_out;
    elist_in_p=&elist_in;
  } else {
    modeout=modein=IGRAPH_ALL;
    IGRAPH_CHECK(igraph_adjedgelist_init(graph,&elist_out, IGRAPH_ALL));
    IGRAPH_FINALLY(igraph_adjedgelist_destroy, &elist_out);
    elist_out_p=elist_in_p=&elist_out;
  }
  
  distance=igraph_Calloc(no_of_nodes, long int);
  if (distance==0) {
    IGRAPH_ERROR("edge betweenness failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, distance);
  nrgeo=igraph_Calloc(no_of_nodes, long int);
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

    IGRAPH_ALLOW_INTERRUPTION();

    memset(distance, 0, no_of_nodes*sizeof(long int));
    memset(nrgeo, 0, no_of_nodes*sizeof(long int));
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
  igraph_vector_t already_counted;
  long int i, j;
  long int nodes_reached;

  igraph_dqueue_t q;
  
  long int nodes_to_calc;
  igraph_vector_t tmp;
  igraph_vit_t vit;

  IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);

  nodes_to_calc=IGRAPH_VIT_SIZE(vit);
  
  if (mode != IGRAPH_OUT && mode != IGRAPH_IN && 
      mode != IGRAPH_ALL) {
    IGRAPH_ERROR("calculating closeness", IGRAPH_EINVMODE);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&already_counted, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);

  IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));
  igraph_vector_null(res);
  
  for (IGRAPH_VIT_RESET(vit), i=0; 
       !IGRAPH_VIT_END(vit); 
       IGRAPH_VIT_NEXT(vit), i++) {
    IGRAPH_CHECK(igraph_dqueue_push(&q, IGRAPH_VIT_GET(vit)));
    IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
    nodes_reached=1;
    VECTOR(already_counted)[(long int)IGRAPH_VIT_GET(vit)]=i+1;

    IGRAPH_ALLOW_INTERRUPTION();
    
    while (!igraph_dqueue_empty(&q)) {
      long int act=igraph_dqueue_pop(&q);
      long int actdist=igraph_dqueue_pop(&q);
      
      VECTOR(*res)[i] += actdist;

      if (cutoff>0 && actdist>=cutoff) continue;

      IGRAPH_CHECK(igraph_neighbors(graph, &tmp, act, mode));
      for (j=0; j<igraph_vector_size(&tmp); j++) {
	long int neighbor=VECTOR(tmp)[j];
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
  
  /* Clean */
  igraph_dqueue_destroy(&q);
  igraph_vector_destroy(&tmp);
  igraph_vector_destroy(&already_counted);
  igraph_vit_destroy(&vit);
  IGRAPH_FINALLY_CLEAN(4);
  
  return 0;
}


