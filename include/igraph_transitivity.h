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

#ifndef IGRAPH_TRANSITIVITY_H
#define IGRAPH_TRANSITIVITY_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

int igraph_transitivity_undirected(const igraph_t *graph, 
				   igraph_real_t *res,
				   igraph_transitivity_mode_t mode);
int igraph_transitivity_local_undirected(const igraph_t *graph, 
					 igraph_vector_t *res,
					 const igraph_vs_t vids,
					 igraph_transitivity_mode_t mode);
int igraph_transitivity_local_undirected1(const igraph_t *graph, 
					  igraph_vector_t *res,
					  const igraph_vs_t vids,
					  igraph_transitivity_mode_t mode);
int igraph_transitivity_local_undirected2(const igraph_t *graph, 
					  igraph_vector_t *res,
					  const igraph_vs_t vids,
					  igraph_transitivity_mode_t mode);
int igraph_transitivity_local_undirected4(const igraph_t *graph, 
					  igraph_vector_t *res,
					  const igraph_vs_t vids,
					  igraph_transitivity_mode_t mode);
int igraph_transitivity_avglocal_undirected(const igraph_t *graph,
					    igraph_real_t *res,
					    igraph_transitivity_mode_t mode);
int igraph_transitivity_barrat(const igraph_t *graph,
			       igraph_vector_t *res,
			       const igraph_vs_t vids,
			       const igraph_vector_t *weights,
			       const igraph_transitivity_mode_t mode);

__END_DECLS

#endif
