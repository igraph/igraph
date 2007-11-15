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

