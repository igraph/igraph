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

int igraph_eigenvector_centrality(const igraph_t *graph, igraph_vector_t *vector,
				  igraph_real_t *value, igraph_bool_t scale,
				  igraph_arpack_options_t *options) {

  igraph_adjlist_t adjlist;
  igraph_vector_t values;
  igraph_matrix_t vectors;
  
  options->n=igraph_vcount(graph);

  IGRAPH_VECTOR_INIT_FINALLY(&values, 0);
  IGRAPH_MATRIX_INIT_FINALLY(&vectors, 0, 0);

  IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

  IGRAPH_CHECK(igraph_arpack_rssolve(igraph_i_eigenvector_centrality,
				     &adjlist, options, &values, &vectors));

  igraph_adjlist_destroy(&adjlist);
  IGRAPH_FINALLY_CLEAN(1);

  if (value) {
    *value=VECTOR(values)[0];
  }
  
  if (vector) {
    igraph_real_t amax=0;
    long int which;
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
				     options, &values, &vectors));

  igraph_adjlist_destroy(&myoutadjlist);
  igraph_adjlist_destroy(&myinadjlist);
  igraph_vector_destroy(&tmp);
  IGRAPH_FINALLY_CLEAN(3);

  if (value) { 
    *value = VECTOR(values)[0];
  }

  if (vector) {
    igraph_real_t amax=0;
    long int which;
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

