/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2013  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_embedding.h"
#include "igraph_interface.h"
#include "igraph_adjlist.h"
#include "igraph_random.h"

typedef struct {
  const igraph_vector_t *cvec;
  igraph_adjlist_t *outlist, *inlist;
  igraph_vector_t *tmp;
} igraph_i_asembedding_data_t;

int igraph_i_asembedding(igraph_real_t *to, const igraph_real_t *from,
			 int n, void *extra) {
  igraph_i_asembedding_data_t *data=extra;
  igraph_adjlist_t *outlist=data->outlist;
  igraph_adjlist_t *inlist=data->inlist;
  const igraph_vector_t *cvec=data->cvec;
  igraph_vector_t *tmp=data->tmp;
  igraph_vector_int_t *neis;
  int i, j, nlen;

  /* tmp = (A+cD)' from */
  for (i=0; i<n; i++) {
    neis=igraph_adjlist_get(inlist, i);
    nlen=igraph_vector_int_size(neis);
    VECTOR(*tmp)[i]=0.0;
    for (j=0; j<nlen; j++) {
      long int nei=(long int) VECTOR(*neis)[j];
      VECTOR(*tmp)[i] += from[nei];
    }
    VECTOR(*tmp)[i] += VECTOR(*cvec)[i] * from[i];
  }
	
  /* to = (A+cD) tmp */
  for (i=0; i<n; i++) {
    neis=igraph_adjlist_get(outlist, i);
    nlen=igraph_vector_int_size(neis);
    to[i]=0.0;
    for (j=0; j<nlen; j++) {
      long int nei=(long int) VECTOR(*neis)[j];
      to[i] += VECTOR(*tmp)[nei];
    }
    to[i] += VECTOR(*cvec)[i] * VECTOR(*tmp)[i];
  }
	
  return 0;
}

int igraph_adjacency_spectral_embedding(const igraph_t *graph, 
					igraph_integer_t no,
					igraph_vector_t *D,
					igraph_matrix_t *U, 
					igraph_matrix_t *V,
					const igraph_vector_t *cvec,
					igraph_arpack_options_t *options) {

  igraph_integer_t vc=igraph_vcount(graph);
  igraph_vector_t tmp;
  igraph_adjlist_t outlist, inlist;
  int i, cveclen=igraph_vector_size(cvec);
  igraph_i_asembedding_data_t data={ cvec, &outlist, &inlist, &tmp };
	
  if (no > vc) { 
    IGRAPH_ERROR("Too many singular values requested", IGRAPH_EINVAL);
  }
  if (no <= 0) {
    IGRAPH_ERROR("No singular values requested", IGRAPH_EINVAL);
  }

  if (cveclen != 1 && cveclen != vc) {
    IGRAPH_ERROR("Augmentation vector size is invalid, it should be "
		 "the number of vertices or scalar", IGRAPH_EINVAL);
  }

  IGRAPH_CHECK(igraph_vector_resize(D, no));
  IGRAPH_CHECK(igraph_matrix_resize(U, vc, no));
  IGRAPH_CHECK(igraph_matrix_resize(V, vc, no));

  /* empty graph */
  if (vc == 1 && igraph_ecount(graph) == 0) {
    igraph_vector_null(D);
    igraph_matrix_null(U);
    igraph_matrix_null(V);
    return 0;
  }

  igraph_vector_init(&tmp, vc);
  IGRAPH_FINALLY(igraph_vector_destroy, &tmp);
  IGRAPH_CHECK(igraph_adjlist_init(graph, &outlist, IGRAPH_OUT));
  IGRAPH_CHECK(igraph_adjlist_init(graph, &inlist, IGRAPH_IN));

  options->n=vc;
  options->start=0;		/* random start vector */
  options->nev=no;
  options->which[0]='L'; options->which[1]='A';
  if (options->ncv <= options->nev) { options->ncv = 0; }
	
  IGRAPH_CHECK(igraph_arpack_rssolve(igraph_i_asembedding,
				     &data, options, 0, D, U));

  if (igraph_is_directed(graph)) {
    data.outlist = &inlist;
    data.inlist  = &outlist;
    IGRAPH_CHECK(igraph_arpack_rssolve(igraph_i_asembedding,
				       &data, options, 0, D, V));
  } else {
    IGRAPH_CHECK(igraph_matrix_update(V, U));
  }
	
  for (i=0; i<no; i++) {
    VECTOR(*D)[i] = sqrt(VECTOR(*D)[i]);
  }
	
  igraph_vector_destroy(&tmp);
  IGRAPH_FINALLY_CLEAN(1);
	
  return 0;
}

igraph_real_t igraph_i_vector_mean(const igraph_vector_t *v, int from,
				   int to) {
  igraph_real_t res=0.0;
  int n=to-from+1;
  for (; from <= to; from++) { res += VECTOR(*v)[from]; }
  return res / n;
}

igraph_real_t igraph_i_vector_var(const igraph_vector_t *v, int from,
				  int to, igraph_real_t mean) {

  igraph_real_t res=0.0;
  int n=to-from+1;
  if (n==1) return 0.0;
  for (; from <= to; from++) {
    res += (mean-VECTOR(*v)[from]) * (mean-VECTOR(*v)[from]);
  }
  return res / (n-1);
}

igraph_real_t igraph_i_vector_log_dnorm(const igraph_vector_t *v, int from,
					int to, igraph_real_t mean,
					igraph_real_t sd) {
  igraph_real_t res=0.0;
  for (; from <= to; from++) {
    res += igraph_dnorm(VECTOR(*v)[from], mean, sd, /*give_log=*/ 1);
  }

  return res;
}

int igraph_dim_select(const igraph_vector_t *sv, igraph_integer_t *dim) {

  int n=igraph_vector_size(sv);
  int i;
  igraph_real_t mean1, mean2, var1, var2, sd, profile;
  igraph_real_t max=IGRAPH_NEGINFINITY;

  if (n==0) {
    IGRAPH_ERROR("Need at least one singular value for dimensionality "
		 "selection", IGRAPH_EINVAL);
  }
  if (n==1) { *dim=1; return 0; }

  for (i = 0; i < n-1; i++) {
    mean1 = igraph_i_vector_mean(sv, 0, i);
    mean2 = igraph_i_vector_mean(sv, i+1, n-1);
    var1 = igraph_i_vector_var(sv, 0, i, mean1);
    var2 = igraph_i_vector_var(sv, i+1, n-1, mean2);
    sd = sqrt((i * var1 + (n-i-2) * var2) / (n-2));
    profile =
      igraph_i_vector_log_dnorm(sv, 0, i, mean1, sd) +
      igraph_i_vector_log_dnorm(sv, i+1, n-1, mean2, sd);
    if (profile > max) {
      max=profile;
      *dim=i+1;
    }
  }

  mean1=igraph_i_vector_mean(sv, 0, n-1);
  sd=sqrt(igraph_i_vector_var(sv, 0, n-1, mean1));
  profile=igraph_i_vector_log_dnorm(sv, 0, n-1, mean1, sd);
  if (profile > max) {
    max=profile;
    *dim=n;
  }

  return 0;
}
