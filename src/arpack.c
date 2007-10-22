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
				  igraph_integer_t maxit, igraph_integer_t pncv) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int i;
  igraph_vector_t *neis;
  igraph_adjlist_t adjlist;

  int n=no_of_nodes, nev=1, ncv=pncv, ldv=n;
  
  igraph_real_t *v, *workl, *workd, *d, *resid, *ax;
  int *select;
  int iparam[11], ipntr[11];
  
  char *bmat="I", *which="LA", *all="All";
  int ido, lworkl, info, ierr, j, ishfts, mode1, nconv, rvec=1;
  igraph_real_t sigma;
  
  igraph_real_t zero=0.0;

  lworkl = ncv*(ncv+8);
  info = 0;
  ido = 0;
  
  ishfts = 1;
  mode1 = 1;
  
  iparam[0]=ishfts;
  iparam[2]=maxit;
  iparam[6]=mode1;

  IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
  v=      igraph_Calloc( ldv*ncv,     igraph_real_t );
  if (!v) {
    IGRAPH_ERROR("Cannot calculate eigenvector centrality", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, v);
  workl=  igraph_Calloc( ncv*(ncv+8), igraph_real_t );
  if (!workl) {
    IGRAPH_ERROR("Cannot calculate eigenvector centrality", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, workl);
  workd=  igraph_Calloc( 3*n,         igraph_real_t );
  if (!workd) {
    IGRAPH_ERROR("Cannot calculate eigenvector centrality", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, workd);
  d=      igraph_Calloc( ncv*2,       igraph_real_t );
  if (!d) {
    IGRAPH_ERROR("Cannot calculate eigenvector centrality", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, d);
  resid=  igraph_Calloc( n,           igraph_real_t );
  if (!resid) {
    IGRAPH_ERROR("Cannot calculate eigenvector centrality", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, resid);
  ax=     igraph_Calloc( n,           igraph_real_t );
  if (!ax) {
    IGRAPH_ERROR("Cannot calculate eigenvector centrality", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, ax);
  select= igraph_Calloc( ncv,         int);
  if (!select) {
    IGRAPH_ERROR("Cannot calculate eigenvector centrality", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, select);

  while (1) {
    
    IGRAPH_F77(dsaupd, DSAUPD) (&ido, bmat, &n, which, &nev, &tol, resid,
				&ncv, v, &ldv, iparam, ipntr, workd, workl,
				&lworkl, &info);
    
    if (ido==-1 || ido==1) {
      
      int j, nlen;
      igraph_real_t *from=&workd[ipntr[0]-1];
      igraph_real_t *to=&workd[ipntr[1]-1];

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

  if (info < 0) {
/*     fprintf(stderr, "ARPACK error %i\n", info); */
    IGRAPH_ERROR("ARPACK error", IGRAPH_FAILURE);
  }

  iparam[4]=1;			/* get eigenvector */
  IGRAPH_F77(dseupd, DSEUPD) (&rvec, all, select, d, v, &ldv, &sigma,
			      bmat, &n, which, &nev, &tol, resid, &ncv,
			      v, &ldv, iparam, ipntr, workd, workl, &lworkl,
			      &ierr);

  if (ierr < 0) {
/*     fprintf(stderr, "ARPACK error %i\n", ierr); */
    IGRAPH_ERROR("ARPACK error", IGRAPH_FAILURE);
  }

  if (value) { 
    *value = d[0];
  }

  if (vector) {
    igraph_real_t amax=0;
    long int which=0;
    IGRAPH_CHECK(igraph_vector_resize(vector, n));
    for (i=0; i<n; i++) {
      igraph_real_t tmp=fabs(v[i]);
      if (tmp > amax) { amax=tmp; which=i; }
      VECTOR(*vector)[i]=v[i];
    }
    if (norm && amax != 0) { igraph_vector_multiply(vector, 1/v[which]); }
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
  igraph_free(ax);
  igraph_free(resid);
  igraph_free(d);
  igraph_free(workd);
  igraph_free(workl);
  igraph_free(v);
  igraph_adjlist_destroy(&adjlist);
  IGRAPH_FINALLY_CLEAN(8);
  
  return 0;
}
