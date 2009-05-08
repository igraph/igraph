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

#ifndef IGRAPH_NONGRAPH_H
#define IGRAPH_NONGRAPH_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

#include "igraph_constants.h"
#include "igraph_types.h"
#include "igraph_matrix.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Other, not graph related                           */
/* -------------------------------------------------- */

int igraph_running_mean(const igraph_vector_t *data, igraph_vector_t *res, 
			igraph_integer_t binwidth);
int igraph_random_sample(igraph_vector_t *res, igraph_integer_t l, igraph_integer_t h, 
			 igraph_integer_t length);
int igraph_convex_hull(const igraph_matrix_t *data, igraph_vector_t *resverts,
		       igraph_matrix_t *rescoords);
int igraph_zeroin(igraph_real_t *ax, igraph_real_t *bx,
		  igraph_real_t (*f)(igraph_real_t x, void *info),
		  void *info, igraph_real_t *Tol, int *Maxit, igraph_real_t *res);
int igraph_bfgs(igraph_vector_t *b, igraph_real_t *Fmin, 
		igraph_scalar_function_t fminfn, igraph_vector_function_t fmingr,
		int maxit, int trace,
		igraph_real_t abstol, igraph_real_t reltol, int nREPORT, void *ex,
		igraph_integer_t *fncount, igraph_integer_t *grcount);

__END_DECLS

#endif
