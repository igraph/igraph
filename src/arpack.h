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

#include "types.h"

#ifndef ARPACK_H
#define ARPACK_H

typedef struct igraph_arpack_options_t {
  /* INPUT */
  char bmat[1];			/* I-standard problem, G-generalized */
  long int n; 			/* Dimension of the eigenproblem */
  char which[2];		/* LA, SA, LM, SM, BE */
  long int nev;                 /* Number of eigenvalues to be computed */
  igraph_real_t tol;		/* Stopping criterion */
  long int ncv;			/* Number of columns in V */
  long int ldv;			/* Leading dimension of V */
  long int ishift;		/* 0-reverse comm., 1-exact with tridiagonal */
  long int mxiter;              /* Maximum number of update iterations to take */
  long int nb;			/* Block size on the recurrence, only 1 works */
  long int mode;		/* The kind of problem to be solved (1-5)
				   1: A*x=l*x, A symmetric
				   2: A*x=l*M*x, A symm. M pos. def.
				   3: K*x = l*M*x, K symm., M pos. semidef.
				   4: K*x = l*KG*x, K s. pos. semidef. KG s. indef.
				   5: A*x = l*M*x, A symm., M symm. pos. semidef. */
  long int start;		/* 0: random, 1: use the supplied vector */
  long int lworkl;		/* Size of temporary storage, default is fine */
  igraph_real_t sigma;          /* The shift for modes 3,4,5 */  
  /* OUTPUT */
  long int info;		/* What happened, see docs */
  long int ierr;		/* What happened  in the dseupd call */
  long int noiter;		/* The number of iterations taken */
  long int numop;		/* Number of OP*x operations */
  long int numopb;		/* Number of B*x operations if BMAT='G' */
  long int numreo;		/* Number of steps of re-orthogonalizations */

  /* INTERNAL */
  long int iparam[11];
  long int ipntr[11];
} igraph_arpack_options_t;

typedef struct igraph_arpack_storage_t {
  long int maxn, maxncv, maxldv;
  igraph_real_t *v;
  igraph_real_t *workl;
  igraph_real_t *workd;
  igraph_real_t *d;
  igraph_real_t *resid;
  igraph_real_t *ax;
  long int *select;
} igraph_arpack_storage_t;

void igraph_arpack_options_init(igraph_arpack_options_t *o);  

int igraph_arpack_storage_init(igraph_arpack_storage_t *s, long int maxn,
			       long int maxncv, long int maxldv);
void igraph_arpack_storage_destroy(igraph_arpack_storage_t *s);

typedef int igraph_arpack_function_t(igraph_real_t *to, const igraph_real_t *from,
				     long int n, void *extra);

int igraph_arpack_rssolve(igraph_arpack_function_t *fun, void *extra,
			  igraph_arpack_options_t *options, 
			  igraph_arpack_storage_t *storage,
			  igraph_vector_t *values, igraph_matrix_t *vectors);

#endif
