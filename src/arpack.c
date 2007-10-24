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
  o->iparam[4]=0; o->iparam[5]=0; o->iparam[6]=o->mode; o->iparam[7]==0;
  o->iparam[8]=0; o->iparam[9]=0; o->iparam[10]=0;
}

int igraph_arpack_rssolve(igraph_arpack_function_t *fun, void *extra,
			  igraph_arpack_options_t *options, 
			  igraph_vector_t *values, igraph_matrix_t *vectors) {
  
  igraph_vector_t v, workl, workd, d, resid, ax; /* just in case */
  int *select;
  igraph_bool_t bv=0, bworkl=0, bworkd=0, bd=0, bresid=0, bax=0, bselect=0;

  int ido=0;
  int rvec= vectors ? 1 : 0;	/* calculate eigenvectors? */
  const char *all="All";
  
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
    select=igraph_Calloc(options->ncv, int);
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
    IGRAPH_F77(dsaupd, DSAUPD)(&ido, options->bmat, &options->n, options->which,
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
    fprintf(stderr, "ARPACK error: %i\n", options->info);
    IGRAPH_ERROR("ARPACK error", IGRAPH_FAILURE);
  }
  
  IGRAPH_F77(dseupd,DSEUPD) (&rvec, all, options->select, VECTOR(*(options->d)),
			     VECTOR(*(options->v)), &options->ldv,
			     &options->sigma, options->bmat, &options->n,
			     options->which, &options->nev, &options->tol,
			     VECTOR(*(options->resid)), &options->ncv,
			     VECTOR(*(options->v)), &options->ldv, options->iparam,
			     options->ipntr, VECTOR(*(options->workd)), 
			     VECTOR(*(options->workl)), &options->lworkl,
			     &options->ierr);
  
  if (options->ierr < 0) {
    fprintf(stderr, "ARPACK error: %i\n", options->ierr);
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

int igraph_i_eigenvector_centrality(igraph_real_t *to, const igraph_real_t *from,
				    int n, void *extra) {
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
				  igraph_real_t *value, igraph_real_t scale,
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
  
  igraph_matrix_destroy(&vectors);
  igraph_vector_destroy(&values);
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}

int igraph_i_kleinberg(const igraph_t *graph, igraph_vector_t *vector,
		       igraph_real_t *value, igraph_integer_t *retcode,
		       igraph_integer_t *vmult, igraph_integer_t *aupdate,
		       igraph_bool_t norm, igraph_real_t tol, 
		       igraph_integer_t maxit, igraph_integer_t pncv,
		       int pwhich, int inout) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int i;
  igraph_vector_t *neis;
  igraph_adjlist_t myinadjlist, myoutadjlist;
  igraph_adjlist_t *inadjlist, *outadjlist;

  int n=no_of_nodes, nev=1, ncv=pncv, ldv=n;
  
  igraph_vector_t v, workl, workd, d, resid, ax, tmp;
  int *select;
  int iparam[11], ipntr[11];
  
  char *bmat="I", *LA="LA", *LM="LM", *all="All", *which;
  int ido, lworkl, info, ierr, j, ishfts, mode1, nconv, rvec=1;
  igraph_real_t sigma;
  
  igraph_real_t zero=0.0;

  if (pwhich==0) { 
    which=LA; 
  } else if (pwhich==1) {
    which=LM; 
  } else {
    IGRAPH_ERROR("Eigenvector centrality failed, `which' is invalid", 
		 IGRAPH_ENOMEM);
  }
  
  if (ncv <= 1) { 
    IGRAPH_ERROR("`ncv'>=3 is required for eigenvector centrality", IGRAPH_EINVAL);
  }

  if (!igraph_is_directed(graph)) {
    IGRAPH_WARNING("Hub scores called on undirected graph");
  }

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

  lworkl = ncv*(ncv+8);
  info = 0;
  ido = 0;
  
  ishfts = 1;
  mode1 = 1;
  
  iparam[0]=ishfts;
  iparam[2]=maxit;
  iparam[6]=mode1;

  IGRAPH_VECTOR_INIT_FINALLY(&v, ldv*ncv);
  IGRAPH_VECTOR_INIT_FINALLY(&workl, lworkl);
  IGRAPH_VECTOR_INIT_FINALLY(&workd, 3*n);
  IGRAPH_VECTOR_INIT_FINALLY(&d, ncv*2);
  IGRAPH_VECTOR_INIT_FINALLY(&resid, n);
  IGRAPH_VECTOR_INIT_FINALLY(&ax, n);
  IGRAPH_VECTOR_INIT_FINALLY(&tmp, n);

  IGRAPH_CHECK(igraph_adjlist_init(graph, &myinadjlist, IGRAPH_IN));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &myinadjlist);
  IGRAPH_CHECK(igraph_adjlist_init(graph, &myoutadjlist, IGRAPH_OUT));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &myoutadjlist);

  select= igraph_Calloc(ncv, int);
  if (!select) {
    IGRAPH_ERROR("Cannot calculate eigenvector centrality", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, select);

  while (1) {
    
    IGRAPH_F77(dsaupd, DSAUPD) (&ido, bmat, &n, which, &nev, &tol, VECTOR(resid),
				&ncv, VECTOR(v), &ldv, iparam, ipntr, 
				VECTOR(workd), VECTOR(workl),
				&lworkl, &info);
    
    if (ido==-1 || ido==1) {
      
      int j, nlen;
      igraph_real_t *from=VECTOR(workd)+ipntr[0]-1;
      igraph_real_t *to=VECTOR(workd)+ipntr[1]-1;

      for (i=0; i<no_of_nodes; i++) {
	neis=igraph_adjlist_get(inadjlist, i);
	nlen=igraph_vector_size(neis);
	VECTOR(tmp)[i]=0.0;
	for (j=0; j<nlen; j++) {
	  long int nei=VECTOR(*neis)[j];
	  VECTOR(tmp)[i] += from[nei];
	}
      }

      for (i=0; i<no_of_nodes; i++) {
	neis=igraph_adjlist_get(outadjlist, i);
	nlen=igraph_vector_size(neis);
	to[i]=0.0;
	for (j=0; j<nlen; j++) {
	  long int nei=VECTOR(*neis)[j];
	  to[i] += VECTOR(tmp)[nei];
	}
      }      
      
    } else {
      break;
    }
    
  }

  igraph_adjlist_destroy(&myoutadjlist);
  igraph_adjlist_destroy(&myinadjlist);
  IGRAPH_FINALLY_CLEAN(2);

  if (info < 0) {
/*     fprintf(stderr, "ARPACK error %i\n", info); */
    IGRAPH_ERROR("ARPACK error", IGRAPH_FAILURE);
  }

  IGRAPH_F77(dseupd, DSEUPD) (&rvec, all, select, VECTOR(d), VECTOR(v), &ldv, &sigma,
			      bmat, &n, which, &nev, &tol, VECTOR(resid), &ncv,
			      VECTOR(v), &ldv, iparam, ipntr, VECTOR(workd), 
			      VECTOR(workl), &lworkl, &ierr);

  if (ierr < 0) {
/*     fprintf(stderr, "ARPACK error %i\n", ierr); */
    IGRAPH_ERROR("ARPACK error", IGRAPH_FAILURE);
  }

  if (value) { 
    *value = VECTOR(d)[0];
  }

  if (vector) {
    igraph_real_t amax=0;
    long int which=0;
    IGRAPH_CHECK(igraph_vector_resize(vector, n));
    for (i=0; i<n; i++) {
      igraph_real_t tmp=fabs(VECTOR(v)[i]);
      if (tmp > amax) { amax=tmp; which=i; }
      VECTOR(*vector)[i]=VECTOR(v)[i];
    }
    if (norm && amax != 0) { igraph_vector_scale(vector, 1/VECTOR(v)[which]); }
  }
  
  if (retcode) {
    *retcode=info;
  }
  if (vmult) {
    *vmult=iparam[8];
  }
  if (aupdate) {
    *aupdate=iparam[2];
  }
  
  /* Free resources */
  igraph_Free(select);
  igraph_vector_destroy(&tmp);
  igraph_vector_destroy(&ax);
  igraph_vector_destroy(&resid);
  igraph_vector_destroy(&d);
  igraph_vector_destroy(&workd);
  igraph_vector_destroy(&workl);
  igraph_vector_destroy(&v);
  IGRAPH_FINALLY_CLEAN(8);
  
  return 0;
}

int igraph_hub_score(const igraph_t *graph, igraph_vector_t *vector,
		     igraph_real_t *value, igraph_integer_t *retcode,
		     igraph_integer_t *vmult, igraph_integer_t *aupdate,
		     igraph_bool_t norm, igraph_real_t tol, 
		     igraph_integer_t maxit, igraph_integer_t pncv,
		     int pwhich) {

  return igraph_i_kleinberg(graph, vector, value, retcode, vmult,
			    aupdate, norm, tol, maxit, pncv, pwhich, 0);
}
			    
int igraph_authority_score(const igraph_t *graph, igraph_vector_t *vector,
			   igraph_real_t *value, igraph_integer_t *retcode,
			   igraph_integer_t *vmult, igraph_integer_t *aupdate,
			   igraph_bool_t norm, igraph_real_t tol, 
			   igraph_integer_t maxit, igraph_integer_t pncv,
			   int pwhich) {

  return igraph_i_kleinberg(graph, vector, value, retcode, vmult,
			    aupdate, norm, tol, maxit, pncv, pwhich, 1);
}

