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

#ifndef IGRAPH_COCITATION_H
#define IGRAPH_COCITATION_H

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

/* -------------------------------------------------- */
/* Cocitation and other similarity measures           */
/* -------------------------------------------------- */

int igraph_cocitation(const igraph_t *graph, igraph_matrix_t *res, 
		      const igraph_vs_t vids);
int igraph_bibcoupling(const igraph_t *graph, igraph_matrix_t *res, 
		       const igraph_vs_t vids);
int igraph_similarity_jaccard(const igraph_t *graph, igraph_matrix_t *res,
              const igraph_vs_t vids, igraph_neimode_t mode,
			  igraph_bool_t loops);
int igraph_similarity_dice(const igraph_t *graph, igraph_matrix_t *res,
              const igraph_vs_t vids, igraph_neimode_t mode,
			  igraph_bool_t loops);
int igraph_similarity_inverse_log_weighted(const igraph_t *graph,
              igraph_matrix_t *res, const igraph_vs_t vids,
              igraph_neimode_t mode);

__END_DECLS

#endif
