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
  o->which[0]='X'; o->which[1]='X';
  o->nev=1;
  o->tol=0;
  o->ncv=3;
  o->ldv=o->n;			/* will be updated to (real) n */
  o->ishift=1;
  o->mxiter=3000;
  o->nb=1;
  o->mode=1;
  o->start=0;
  o->lworkl=0;
  o->sigma=0;
  o->sigmai=0;
  o->info=o->start;
  
  o->iparam[0]=o->ishift; o->iparam[1]=0; o->iparam[2]=o->mxiter; o->iparam[3]=o->nb;
  o->iparam[4]=0; o->iparam[5]=0; o->iparam[6]=o->mode; o->iparam[7]=0;
  o->iparam[8]=0; o->iparam[9]=0; o->iparam[10]=0;
}

int igraph_arpack_storage_init(igraph_arpack_storage_t *s, long int maxn,
			       long int maxncv, long int maxldv, 
			       igraph_bool_t symm) {
  
  /* TODO: check arguments */
  s->maxn=maxn; 
  s->maxncv=maxncv;
  s->maxldv=maxldv;

#define CHECKMEM(x) \
    if (!x) { \
      IGRAPH_ERROR("Cannot allocate memory for ARPACK", IGRAPH_ENOMEM); \
    } \
    IGRAPH_FINALLY(igraph_free, x);

  s->v=igraph_Calloc(maxldv * maxncv, igraph_real_t); CHECKMEM(s->v);
  s->workd=igraph_Calloc(3*maxn, igraph_real_t); CHECKMEM(s->workd);
  s->d=igraph_Calloc(2*maxncv, igraph_real_t); CHECKMEM(s->d);
  s->resid=igraph_Calloc(maxn, igraph_real_t); CHECKMEM(s->resid);
  s->ax=igraph_Calloc(maxn, igraph_real_t); CHECKMEM(s->ax);
  s->select=igraph_Calloc(maxncv, long int); CHECKMEM(s->select);

  if (symm) {
    s->workl=igraph_Calloc(maxncv*(maxncv+8), igraph_real_t); CHECKMEM(s->workl);
    s->di=0;
    s->workev=0;
  } else {
    s->workl=igraph_Calloc(3*maxncv*(maxncv+2), igraph_real_t); CHECKMEM(s->workl);
    s->di=igraph_Calloc(2*maxncv, igraph_real_t); CHECKMEM(s->di);
    s->workev=igraph_Calloc(3*maxncv, igraph_real_t); CHECKMEM(s->workev);
    IGRAPH_FINALLY_CLEAN(2);
  }

#undef CHECKMEM  
  
  IGRAPH_FINALLY_CLEAN(7);
  return 0;
}

void igraph_arpack_storage_destroy(igraph_arpack_storage_t *s) {
  
  if (s->di) { 
    igraph_Free(s->di);
  }
  if (s->workev) {
    igraph_Free(s->workev);
  }
  
  igraph_Free(s->workl);
  igraph_Free(s->select);
  igraph_Free(s->ax);
  igraph_Free(s->resid);
  igraph_Free(s->d);
  igraph_Free(s->workd);
  igraph_Free(s->v);
}

int igraph_arpack_rssolve(igraph_arpack_function_t *fun, void *extra,
			  igraph_arpack_options_t *options, 
			  igraph_arpack_storage_t *storage,
			  igraph_vector_t *values, igraph_matrix_t *vectors) {
  
  igraph_real_t *v, *workl, *workd, *d, *resid, *ax;
  igraph_bool_t free_them=0;
  long int *select;

  long int ido=0;
  long int rvec= vectors || storage ? 1 : 0;	/* calculate eigenvectors? */
  char *all="All";
  
  /* Brush up options if needed */
  if (options->ldv == 0) { options->ldv=options->n; }
  if (options->lworkl == 0) { options->lworkl=options->ncv*(options->ncv+8); }
  if (options->which[0] == 'X') { options->which[0]='L'; options->which[1]='A'; }
  
  if (storage) {
    /* Storage provided */
    if (storage->maxn < options->n) {
      IGRAPH_ERROR("Not enough storage for ARPACK (`n')", IGRAPH_EINVAL);
    }
    if (storage->maxncv < options->ncv) {
      IGRAPH_ERROR("Not enough storage for ARPACK (`ncv')", IGRAPH_EINVAL);
    }
    if (storage->maxldv < options->ldv) {
      IGRAPH_ERROR("Not enough storage for ARPACK (`ldv')", IGRAPH_EINVAL);
    }

    v      = storage->v; 
    workl  = storage->workl;
    workd  = storage->workd;
    d      = storage->d;
    resid  = storage->resid;
    ax     = storage->ax;
    select = storage->select;
    
  } else {
    /* Storage not provided */
    free_them=1;
    
#define CHECKMEM(x) \
    if (!x) { \
      IGRAPH_ERROR("Cannot allocate memory for ARPACK", IGRAPH_ENOMEM); \
    } \
    IGRAPH_FINALLY(igraph_free, x);

    v=igraph_Calloc(options->ldv * options->ncv, igraph_real_t); CHECKMEM(v);
    workl=igraph_Calloc(options->lworkl, igraph_real_t); CHECKMEM(workl);
    workd=igraph_Calloc(3*options->n, igraph_real_t); CHECKMEM(workd);
    d=igraph_Calloc(2*options->ncv, igraph_real_t); CHECKMEM(d);
    resid=igraph_Calloc(options->n, igraph_real_t); CHECKMEM(resid);
    ax=igraph_Calloc(options->n, igraph_real_t); CHECKMEM(ax);
    select=igraph_Calloc(options->ncv, long int); CHECKMEM(select);

#undef CHECKMEM

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
		  resid, &options->ncv, v, &options->ldv,
		  options->iparam, options->ipntr,
		  workd, workl, &options->lworkl, &options->info);

    if (ido==-1 || ido==1) {

      igraph_real_t *from=workd+options->ipntr[0]-1;
      igraph_real_t *to=workd+options->ipntr[1]-1;
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
  
  igraphdseupd_(&rvec, all, select, d, v, &options->ldv,
		&options->sigma, options->bmat, &options->n,
		options->which, &options->nev, &options->tol,
		resid, &options->ncv, v, &options->ldv, options->iparam,
		options->ipntr, workd, workl, &options->lworkl,
		&options->ierr);
  
  if (options->ierr < 0) {
/*     fprintf(stderr, "ARPACK error: %i\n", (int)options->ierr); */
    IGRAPH_ERROR("ARPACK error", IGRAPH_FAILURE);
  }    
  
  /* Save the result */
  
  if (values) {
    long int i;
    IGRAPH_CHECK(igraph_vector_resize(values, options->nev));
    for (i=0; i<options->nev; i++) {
      VECTOR(*values)[i] = d[i];
    }
  }

  if (vectors) {
    long int i, j, ptr=0;
    IGRAPH_CHECK(igraph_matrix_resize(vectors, options->n, options->nev));
    for (j=0; j<options->nev; j++) {	
      for (i=0; i<options->n; i++) {
	MATRIX(*vectors, i, j) = v[ptr++];
      }
    }
  }
  
  /* Clean up if needed */
  if (free_them) {
    igraph_Free(select);
    igraph_Free(ax);
    igraph_Free(resid);
    igraph_Free(d);
    igraph_Free(workd);
    igraph_Free(workl);
    igraph_Free(v);
    IGRAPH_FINALLY_CLEAN(7);
  }
  return 0;
}

int igraph_arpack_rnsolve(igraph_arpack_function_t *fun, void *extra,
			  igraph_arpack_options_t *options,
			  igraph_arpack_storage_t *storage,
			  igraph_matrix_t *values, igraph_matrix_t *vectors) {

  igraph_real_t *v, *workl, *workd, *dr, *di, *resid, *ax, *workev;
  igraph_bool_t free_them=0;
  long int *select;
  
  long int ido=0;
  long int rvec= vectors ? 1 : 0;
  char *all="All";
  
  /* Brush up options if needed */
  if (options->ldv == 0) { options->ldv=options->n; }
  if (options->lworkl == 0) { options->lworkl=3*options->ncv*(options->ncv+2); }
  if (options->which[0] == 'X') { options->which[0]='L'; options->which[1]='M'; }
  
  if (storage) {
    /* Storage provided */
    if (storage->maxn < options->n) {
      IGRAPH_ERROR("Not enough storage for ARPACK (`n')", IGRAPH_EINVAL);
    }
    if (storage->maxncv < options->ncv) {
      IGRAPH_ERROR("Not enough storage for ARPACK (`ncv')", IGRAPH_EINVAL);
    }
    if (storage->maxldv < options->ldv) {
      IGRAPH_ERROR("Not enough storage for ARPACK (`ldv')", IGRAPH_EINVAL);
    }

    v      = storage->v; 
    workl  = storage->workl;
    workd  = storage->workd;
    workev = storage->workev;
    dr     = storage->d;
    di     = storage->di;
    resid  = storage->resid;
    ax     = storage->ax;
    select = storage->select;
    
  } else {
    /* Storage not provided */
    free_them=1;
    
#define CHECKMEM(x) \
    if (!x) { \
      IGRAPH_ERROR("Cannot allocate memory for ARPACK", IGRAPH_ENOMEM); \
    } \
    IGRAPH_FINALLY(igraph_free, x);

    v=igraph_Calloc(options->ldv * options->ncv, igraph_real_t); CHECKMEM(v);
    workl=igraph_Calloc(options->lworkl, igraph_real_t); CHECKMEM(workl);
    workd=igraph_Calloc(3*options->n, igraph_real_t); CHECKMEM(workd);
    dr=igraph_Calloc(2*options->ncv, igraph_real_t); CHECKMEM(dr);
    di=igraph_Calloc(2*options->ncv, igraph_real_t); CHECKMEM(di);
    resid=igraph_Calloc(options->n, igraph_real_t); CHECKMEM(resid);
    ax=igraph_Calloc(options->n, igraph_real_t); CHECKMEM(ax);
    select=igraph_Calloc(options->ncv, long int); CHECKMEM(select);
    workev=igraph_Calloc(3*options->ncv, igraph_real_t); CHECKMEM(workev);
    
#undef CHECKMEM

  }

  /* Set final bits */
  options->iparam[0]=options->ishift;
  options->iparam[2]=options->mxiter;
  options->iparam[6]=options->mode;
  options->info=options->start;
  
  /* Ok, we have everything */
  while (1) {
    
    igraphdnaupd_(&ido, options->bmat, &options->n, options->which,
		  &options->nev, &options->tol,
		  resid, &options->ncv, v, &options->ldv,
		  options->iparam, options->ipntr,
		  workd, workl, &options->lworkl, &options->info);
    
    if (ido==-1 || ido==1) {
      igraph_real_t *from=workd+options->ipntr[0]-1;
      igraph_real_t *to=workd+options->ipntr[1]-1;
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

  igraphdneupd_(&rvec, all, select, dr, di, v, &options->ldv,
		&options->sigma, &options->sigmai, workev, options->bmat,
		&options->n, options->which, &options->nev, &options->tol,
		resid, &options->ncv, v, &options->ldv, options->iparam,
		options->ipntr, workd, workl, &options->lworkl,
		&options->ierr);
  
  if (options->ierr < 0) {
/*     fprintf(stderr, "ARPACK error: %i\n", (int)options->ierr); */
    IGRAPH_ERROR("ARPACK error", IGRAPH_FAILURE);
  }    

  /* Save the result */
  
  if (values) {
    long int i;
    IGRAPH_CHECK(igraph_matrix_resize(values, options->nev, 2));
    for (i=0; i<options->nev; i++) {
      MATRIX(*values, i, 0) = dr[i];
      MATRIX(*values, i, 1) = di[i];
    }
  }
  
  if (vectors) {
    long int i, j, ptr=0;
    IGRAPH_CHECK(igraph_matrix_resize(vectors, options->n, options->nev+1));
    for (j=0; j<options->nev+1; j++) {	
      for (i=0; i<options->n; i++) {
	MATRIX(*vectors, i, j) = v[ptr++];
      }
    }
  }
  
  /* Clean up if needed */
  if (free_them) {
    igraph_Free(workev);
    igraph_Free(select);
    igraph_Free(ax);
    igraph_Free(resid);
    igraph_Free(di);
    igraph_Free(dr);
    igraph_Free(workd);
    igraph_Free(workl);
    igraph_Free(v);
    IGRAPH_FINALLY_CLEAN(8);
  }
  return 0;  
}
