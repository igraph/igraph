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

#include "arpack.h"
#include "arpack_internal.h"
#include "memory.h"
#include "igraph.h"

#include <math.h>
#include <stdio.h>

/* The ARPACK example file dssimp.f is used as a template */

void igraph_arpack_options_init(igraph_arpack_options_t *o) {
  o->bmat[0]='I';
  o->n=0;			/* needs to be updated! */
  o->which[0]='L'; o->which[1]='A';
  o->nev=1;
  o->tol=0;
  o->ncv=3;
  o->ldv=o->n;			/* will be updated to (real) n */
  o->ishift=1;
  o->mxiter=300;
  o->nb=1;
  o->mode=1;
  o->start=0;
  o->lworkl=o->ncv*(o->ncv+8);
  o->sigma=0;
  o->info=o->start;
  
  o->v=0;
  o->workl=0;
  o->workd=0;
  o->d=0;
  o->resid=0;
  o->ax=0;
  o->select=0;
  o->iparam[0]=o->ishift; o->iparam[1]=0; o->iparam[2]=o->mxiter; o->iparam[3]=o->nb;
  o->iparam[4]=0; o->iparam[5]=0; o->iparam[6]=o->mode; o->iparam[7]=0;
  o->iparam[8]=0; o->iparam[9]=0; o->iparam[10]=0;
}

int igraph_arpack_rssolve(igraph_arpack_function_t *fun, void *extra,
			  igraph_arpack_options_t *options, 
			  igraph_vector_t *values, igraph_matrix_t *vectors) {
  
  igraph_vector_t v, workl, workd, d, resid, ax; /* just in case */
  long int *select=0;
  igraph_bool_t bv=0, bworkl=0, bworkd=0, bd=0, bresid=0, bax=0, bselect=0;

  long int ido=0;
  long int rvec= vectors ? 1 : 0;	/* calculate eigenvectors? */
  char *all="All";
  
  /* Brush up options if needed */
  if (options->ldv == 0) { options->ldv=options->n; }
  if (options->lworkl == 0) { options->lworkl=options->ncv*(options->ncv+8); }
  if (!options->v) {
    bv=1;
    IGRAPH_VECTOR_INIT_FINALLY(&v, options->ldv * options->ncv);
    options->v=&v;
  }
  if (!options->workl) {
    bworkl=1;
    IGRAPH_VECTOR_INIT_FINALLY(&workl, options->lworkl);
    options->workl=&workl;
  }
  if (!options->workd) {
    bworkd=1;
    IGRAPH_VECTOR_INIT_FINALLY(&workd, 3*options->n);
    options->workd=&workd;
  }
  if (!options->d) {
    bd=1;
    IGRAPH_VECTOR_INIT_FINALLY(&d, 2*options->ncv);
    options->d=&d;
  }
  if (!options->resid) {
    bresid=1;
    IGRAPH_VECTOR_INIT_FINALLY(&resid, options->n);
    options->resid=&resid;
  }
  if (!options->ax) {
    bax=1;
    IGRAPH_VECTOR_INIT_FINALLY(&ax, options->n);
    options->ax=&ax;
  }
  if (!options->select) {
    bselect=1;
    select=igraph_Calloc(options->ncv, long int);
    if (!select) { 
      IGRAPH_ERROR("Cannot do rssolve", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, select);
    options->select=select;
  }

  /* Set final bits */
  options->iparam[0]=options->ishift;
  options->iparam[2]=options->mxiter;
  options->iparam[6]=options->mode;
  options->info=options->start;
  
  /* Ok, we have everything */
  while (1) {
    igraphdsaupd_(&ido, options->bmat, &options->n, options->which,
		  &options->nev, &options->tol, 
		  VECTOR(*(options->resid)), &options->ncv,
		  VECTOR(*(options->v)), &options->ldv, 
		  options->iparam, options->ipntr,
		  VECTOR(*(options->workd)), VECTOR(*(options->workl)),
		  &options->lworkl, &options->info);
    if (ido==-1 || ido==1) {

      igraph_real_t *from=VECTOR(*(options->workd))+options->ipntr[0]-1;
      igraph_real_t *to=VECTOR(*(options->workd))+options->ipntr[1]-1;
      if (fun(to, from, options->n, extra) != 0) {
	IGRAPH_ERROR("Arpack error while evaluating matrix-vector product",
		     IGRAPH_FAILURE);
      }
      
    } else {
      break;
    }
  }
  
  if (options->info < 0) {
    fprintf(stderr, "ARPACK error: %i\n", (int) options->info);
    IGRAPH_ERROR("ARPACK error", IGRAPH_FAILURE);
  }
  
  igraphdseupd_(&rvec, all, options->select, VECTOR(*(options->d)),
		VECTOR(*(options->v)), &options->ldv,
		&options->sigma, options->bmat, &options->n,
		options->which, &options->nev, &options->tol,
		VECTOR(*(options->resid)), &options->ncv,
		VECTOR(*(options->v)), &options->ldv, options->iparam,
		options->ipntr, VECTOR(*(options->workd)), 
		VECTOR(*(options->workl)), &options->lworkl,
		&options->ierr);
  
  if (options->ierr < 0) {
    fprintf(stderr, "ARPACK error: %i\n", (int)options->ierr);
    IGRAPH_ERROR("ARPACK error", IGRAPH_FAILURE);
  }    

  /* Save the result */
  
  if (values) {
    long int i;
    IGRAPH_CHECK(igraph_vector_resize(values, options->nev));
    for (i=0; i<options->nev; i++) {
      VECTOR(*values)[i] = VECTOR(*(options->d))[i];
    }
  }
  
  if (vectors) {
    long int i, j, ptr=0;
    IGRAPH_CHECK(igraph_matrix_resize(vectors, options->n, options->nev));
    for (i=0; i<options->n; i++) {
      for (j=0; j<options->nev; j++) {	
	MATRIX(*vectors, i, j) = VECTOR(*(options->v))[ptr++];
      }
    }
  }
  
  /* Clean up */
  if (bselect) {
    igraph_Free(select);
    options->select=0;
    IGRAPH_FINALLY_CLEAN(1);
  }
  if (bax) {
    igraph_vector_destroy(&ax);
    options->ax=0;
    IGRAPH_FINALLY_CLEAN(1);
  }
  if (bresid) {
    igraph_vector_destroy(&resid);
    options->resid=0;
    IGRAPH_FINALLY_CLEAN(1);
  }
  if (bd) {
    igraph_vector_destroy(&d);
    options->d=0;
    IGRAPH_FINALLY_CLEAN(1);
  }
  if (bworkd) {
    igraph_vector_destroy(&workd);
    options->workd=0;
    IGRAPH_FINALLY_CLEAN(1);
  }
  if (bworkl) {
    igraph_vector_destroy(&workl);
    options->workl=0;
    IGRAPH_FINALLY_CLEAN(1);
  }
  return 0;
}
