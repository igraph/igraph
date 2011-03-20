/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2009  Gabor Csardi <csardi@rmki.kfki.hu>
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

#ifndef IGRAPH_HRG_H
#define IGRAPH_HRG_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

#include "igraph_vector.h"
#include "igraph_vector_ptr.h"
#include "igraph_datatype.h"

__BEGIN_DECLS

typedef struct igraph_hrg_t {
  igraph_vector_t left, right, prob, edges, vertices;
} igraph_hrg_t;

int igraph_hrg_init(igraph_hrg_t *hrg, int no_of_nodes);
void igraph_hrg_destroy(igraph_hrg_t *hrg);
int igraph_hrg_size(const igraph_hrg_t *hrg);
int igraph_hrg_resize(igraph_hrg_t *hrg, int newsize);

int igraph_hrg_fit(const igraph_t *graph, 
		   igraph_hrg_t *hrg,
		   igraph_bool_t start,
		   int steps);

int igraph_hrg_sample(const igraph_t *graph,
		      igraph_t *sample,
		      igraph_vector_ptr_t *samples,
		      igraph_hrg_t *hrg,
		      igraph_bool_t start);

int igraph_hrg_game(igraph_t *graph,
		    const igraph_hrg_t *hrg);

int igraph_hrg_dendrogram(igraph_t *graph,
			  const igraph_hrg_t *hrg);

int igraph_hrg_consensus(const igraph_t *graph,
			 igraph_vector_t *parents,
			 igraph_vector_t *weights,
			 igraph_hrg_t *hrg,
			 igraph_bool_t start, 
			 int num_samples);

int igraph_hrg_predict(const igraph_t *graph,
		       igraph_vector_t *edges,
		       igraph_vector_t *prob,
		       igraph_hrg_t *hrg,
		       igraph_bool_t start, 
		       int num_samples,
		       int num_bins);

int igraph_hrg_create(igraph_hrg_t *hrg,
		      const igraph_t *graph, 
		      const igraph_vector_t *prob);

__END_DECLS

#endif	/* IGRAPH_HRG_H */
