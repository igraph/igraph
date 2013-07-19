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

typedef struct {
	igraph_real_t c;
	igraph_vector_t *degree;
	igraph_adjlist_t *outlist, *inlist;
	igraph_vector_t *tmp;
} igraph_i_asembedding_data_t;

int igraph_i_asembedding(igraph_real_t *to, const igraph_real_t *from,
												 int n, void *extra) {
	igraph_i_asembedding_data_t *data=extra;
	igraph_real_t c=data->c;
	igraph_adjlist_t *outlist=data->outlist;
	igraph_adjlist_t *inlist=data->inlist;
	igraph_vector_t *degree=data->degree;
	igraph_vector_t *tmp=data->tmp;
	igraph_vector_t *neis;
	int i, j, nlen;

	/* tmp = (A+cD)' from */
	for (i=0; i<n; i++) {
		neis=igraph_adjlist_get(inlist, i);
		nlen=igraph_vector_size(neis);
		VECTOR(*tmp)[i]=0.0;
		for (j=0; j<nlen; j++) {
			long int nei=(long int) VECTOR(*neis)[j];
			VECTOR(*tmp)[i] += from[nei];
		}
		VECTOR(*tmp)[i] += c * VECTOR(*degree)[i] * from[i];
	}
	
	/* to = (A+cD) tmp */
	for (i=0; i<n; i++) {
		neis=igraph_adjlist_get(outlist, i);
		nlen=igraph_vector_size(neis);
		to[i]=0.0;
		for (j=0; j<nlen; j++) {
			long int nei=(long int) VECTOR(*neis)[j];
			to[i] += VECTOR(*tmp)[nei];
		}
		to[i] += c * VECTOR(*degree)[i] * VECTOR(*tmp)[i];
	}
	
	return 0;
}

int igraph_adjacency_spectral_embedding(const igraph_t *graph, 
																				igraph_integer_t no,
																				igraph_vector_t *D,
																				igraph_matrix_t *U, 
																				igraph_matrix_t *V, 
																				igraph_real_t c,
																				igraph_arpack_options_t *options) {

	igraph_integer_t vc=igraph_vcount(graph);
	igraph_vector_t deg, tmp;
	igraph_adjlist_t outlist, inlist;
	int i;
	igraph_i_asembedding_data_t data={ c, &deg, &outlist, &inlist, &tmp };
	
	if (no > vc) { 
		IGRAPH_ERROR("Too many singular values requested", IGRAPH_EINVAL);
	}
	if (no <= 0) {
		IGRAPH_ERROR("No singular values requested", IGRAPH_EINVAL);
	}

	IGRAPH_CHECK(igraph_vector_resize(D, no));
	IGRAPH_CHECK(igraph_matrix_resize(U, vc, no));
	IGRAPH_CHECK(igraph_matrix_resize(V, vc, no));

	/* empty graph */
	if (igraph_ecount(graph) == 0) {
		igraph_vector_null(D);
		igraph_matrix_null(U);
		igraph_matrix_null(V);
		return 0;
	}

	igraph_vector_init(&tmp, vc);
	IGRAPH_FINALLY(igraph_vector_destroy, &tmp);
	igraph_vector_init(&deg, vc);
	IGRAPH_FINALLY(igraph_vector_destroy, &deg);
	IGRAPH_CHECK(igraph_degree(graph, &deg, igraph_vss_all(), IGRAPH_ALL, 
														 /*loops=*/ 0));
	
	IGRAPH_CHECK(igraph_adjlist_init(graph, &outlist, IGRAPH_OUT));
	IGRAPH_CHECK(igraph_adjlist_init(graph, &inlist, IGRAPH_IN));

	options->n=vc;
	options->start=0;							/* random start vector */
	options->nev=no;
  options->which[0]='L'; options->which[1]='A';
	
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
	
	igraph_vector_destroy(&deg);
	igraph_vector_destroy(&tmp);
	IGRAPH_FINALLY_CLEAN(2);
	
	return 0;
}
