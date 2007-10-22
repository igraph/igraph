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

#include "igraph.h"
#include "config.h"
#define IGRAPH_NO_CALLOC
#include "memory.h"

#include <math.h>

#ifdef USING_R
/* #  include <R.h> */
/* #  define IGRAPH_F77(a,b) F77_CALL(a) */
#  define IGRAPH_F77(a,b) igraph ## a ## _
#else
#  define IGRAPH_F77(a,b) F77_FUNC(igraph ## a, IGRAPH ## b)
#endif

void IGRAPH_F77(dsaupd,DSAUPD)(int *ido, const char *bmat, int *n,
			       const char *which, int *nev, igraph_real_t *tol,
			       igraph_real_t *resid, int *ncv, igraph_real_t *v,
			       int *ldv, int *iparam, int *ipntr, 
			       igraph_real_t *workd, igraph_real_t *workl,
			       int *lworkl, int *info);

void IGRAPH_F77(dseupd,DSEUPD)(int *rvec, const char *howmny, int *select,
			       igraph_real_t *d, igraph_real_t *z, int *ldz,
			       igraph_real_t *sigma, const char *bmat, int *n,
			       const char *which, int *nev, igraph_real_t *tol,
			       igraph_real_t *resid, int *ncv, igraph_real_t *v,
			       int *ldv, int *iparam, int *ipntr, 
			       igraph_real_t *workd, igraph_real_t *workl,
			       int *lworkl, int *info);

/* The ARPACK example file dssimp.f is used as a template */

int igraph_eigenvector_centrality(const igraph_t *graph, igraph_vector_t *vector,
				  igraph_real_t *value, igraph_integer_t *retcode,
				  igraph_integer_t *vmult, igraph_integer_t *aupdate,
				  igraph_bool_t norm, igraph_real_t tol, 
				  igraph_integer_t maxit, igraph_integer_t pncv,
				  int pwhich) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int i;
  igraph_vector_t *neis;
  igraph_adjlist_t adjlist;

  int n=no_of_nodes, nev=1, ncv=pncv, ldv=n;
  
  igraph_vector_t v, workl, workd, d, resid, ax;
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
  select= igraph_Calloc(ncv, int);
  if (!select) {
    IGRAPH_ERROR("Cannot calculate eigenvector centrality", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, select);

  IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

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
	neis=igraph_adjlist_get(&adjlist, i);
	nlen=igraph_vector_size(neis);
	to[i]=0.0;
	for (j=0; j<nlen; j++) {
	  long int nei=VECTOR(*neis)[j];
	  to[i] += from[nei];
	}
      }				      
      
    } else {
      break;
    }
    
  }

  igraph_adjlist_destroy(&adjlist);
  IGRAPH_FINALLY_CLEAN(1);

  if (info < 0) {
/*     fprintf(stderr, "ARPACK error %i\n", info); */
    IGRAPH_ERROR("ARPACK error", IGRAPH_FAILURE);
  }

  iparam[4]=1;			/* get eigenvector */
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
    if (norm && amax != 0) { igraph_vector_multiply(vector, 1/VECTOR(v)[which]); }
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
  igraph_free(select);
  igraph_vector_destroy(&ax);
  igraph_vector_destroy(&resid);
  igraph_vector_destroy(&d);
  igraph_vector_destroy(&workd);
  igraph_vector_destroy(&workl);
  igraph_vector_destroy(&v);
  IGRAPH_FINALLY_CLEAN(7);
  
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

  iparam[4]=1;			/* get eigenvector */
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
    if (norm && amax != 0) { igraph_vector_multiply(vector, 1/VECTOR(v)[which]); }
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
  igraph_free(select);
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

