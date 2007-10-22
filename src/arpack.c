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
#include "memory.h"

#define IGRAPH_F77(a,b) F77_FUNC(igraph ## a, IGRAPH ## b)

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

int igraph_eigenvector_centrality(const igraph_t *graph, igraph_vector_t *res) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int i;
  igraph_vector_t neis;

  /* %------------------------------------------------------%
     | Storage Declarations:                                |
     |                                                      |
     | The maximum dimensions for all arrays are            |
     | set here to accommodate a problem size of            |
     | N .le. MAXN                                          |
     |                                                      |
     | NEV is the number of eigenvalues requested.          |
     |     See specifications for ARPACK usage below.       |
     |                                                      |
     | NCV is the largest number of basis vectors that will |
     |     be used in the Implicitly Restarted Arnoldi      |
     |     Process.  Work per major iteration is            |
     |     proportional to N*NCV*NCV.                       |
     |                                                      |
     | You must set:                                        |
     |                                                      |
     | MAXN:   Maximum dimension of the A allowed.          |
     | MAXNEV: Maximum NEV allowed.                         |
     | MAXNCV: Maximum NCV allowed.                         |
     %------------------------------------------------------% */
  
  int maxn=no_of_nodes, maxnev=1, maxncv=3, ldv=maxn;
  
  /* %--------------%
     | Local Arrays |
     %--------------% */

  igraph_real_t *v, *workl, *workd, *d, *resid, *ax;
  int *select;
  int iparam[11], ipntr[11];
  
  /* %---------------%
     | Local Scalars |
     %---------------% */

  char *bmat="I", *which="LM", *all="All";
  int ido, n, nev, ncv, lworkl, info, ierr, j, nx, ishfts, maxitr, 
    mode1, nconv, rvec=1;
  igraph_real_t tol, sigma;
  
  /* %------------%
     | Parameters |
     %------------% */

  igraph_real_t zero=0.0;
  
  /* %-----------------------%
     | Executable Statements |
     %-----------------------% */
  
  /* %-------------------------------------------------%
     | The following include statement and assignments |
     | initiate trace output from the internal         |
     | actions of ARPACK.  See debug.doc in the        |
     | DOCUMENTS directory for usage.  Initially, the  |
     | most useful information will be a breakdown of  |
     | time spent in the various stages of computation |
     | given by setting msaupd = 1.                    |
     %-------------------------------------------------% */

  int  logfil, ndigit, mgetv0,
    msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
    mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
    mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd;

  ndigit = -3;
  logfil = 6;
  msgets = 0;
  msaitr = 0;
  msapps = 0;
  msaupd = 1;
  msaup2 = 0;
  mseigt = 0;
  mseupd = 0;
  
  /* %-------------------------------------------------%
     | The following sets dimensions for this problem. |
     %-------------------------------------------------% */
  
  nx = no_of_nodes;
  n = nx;

  /* %-----------------------------------------------%
     |                                               | 
     | Specifications for ARPACK usage are set       | 
     | below:                                        |
     |                                               |
     |    1) NEV = 4  asks for 4 eigenvalues to be   |  
     |       computed.                               | 
     |                                               |
     |    2) NCV = 20 sets the length of the Arnoldi |
     |       factorization                           |
     |                                               |
     |    3) This is a standard problem              |
     |         (indicated by bmat  = 'I')            |
     |                                               |
     |    4) Ask for the NEV eigenvalues of          |
     |       largest magnitude                       |
     |         (indicated by which = 'LM')           |
     |       See documentation in DSAUPD for the     |
     |       other options SM, LA, SA, LI, SI.       | 
     |                                               |
     | Note: NEV and NCV must satisfy the following  |
     | conditions:                                   |
     |              NEV <= MAXNEV                    |
     |          NEV + 1 <= NCV <= MAXNCV             |
     %-----------------------------------------------% */

  nev=1;
  ncv=3;
  /* bmat is set to "I" */
  /* which is set to "LA" */
  
  /* %-----------------------------------------------------%
     |                                                     |
     | Specification of stopping rules and initial         |
     | conditions before calling DSAUPD                    |
     |                                                     |
     | TOL  determines the stopping criterion.             |
     |                                                     |
     |      Expect                                         |
     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) |
     |               computed   true                       |
     |                                                     |
     |      If TOL .le. 0,  then TOL <- macheps            |
     |           (machine precision) is used.              |
     |                                                     |
     | IDO  is the REVERSE COMMUNICATION parameter         |
     |      used to specify actions to be taken on return  |
     |      from DSAUPD. (See usage below.)                |
     |                                                     |
     |      It MUST initially be set to 0 before the first |
     |      call to DSAUPD.                                | 
     |                                                     |
     | INFO on entry specifies starting vector information |
     |      and on return indicates error codes            |
     |                                                     |
     |      Initially, setting INFO=0 indicates that a     | 
     |      random starting vector is requested to         |
     |      start the ARNOLDI iteration.  Setting INFO to  |
     |      a nonzero value on the initial call is used    |
     |      if you want to specify your own starting       |
     |      vector (This vector must be placed in RESID.)  | 
     |                                                     |
     | The work array WORKL is used in DSAUPD as           | 
     | workspace.  Its dimension LWORKL is set as          |
     | illustrated below.                                  |
     |                                                     |
     %-----------------------------------------------------% */
  
  lworkl = ncv*(ncv+8);
  tol = 0;
  info = 0;
  ido = 0;
  
  /* %---------------------------------------------------%
     | Specification of Algorithm Mode:                  |
     |                                                   |
     | This program uses the exact shift strategy        |
     | (indicated by setting PARAM(1) = 1).              |
     | IPARAM(3) specifies the maximum number of Arnoldi |
     | iterations allowed.  Mode 1 of DSAUPD is used     |
     | (IPARAM(7) = 1). All these options can be changed |
     | by the user. For details see the documentation in |
     | DSAUPD.                                           |
     %---------------------------------------------------% */

  ishfts = 1;
  maxitr = 300;
  mode1 = 1;
  
  iparam[0]=ishfts;
  iparam[2]=maxitr;
  iparam[6]=mode1;

  /* This is not part of the template, we need to allocate memory 
     for the arrays */

  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  v=      Calloc( ldv*maxncv,        igraph_real_t );
  workl=  Calloc( maxncv*(maxncv+8), igraph_real_t );
  workd=  Calloc( 3*maxn,            igraph_real_t );
  d=      Calloc( maxncv*2,          igraph_real_t );
  resid=  Calloc( maxn,              igraph_real_t );
  ax=     Calloc( maxn,              igraph_real_t );
  select= Calloc( maxncv,            int);

  /* %------------------------------------------------%
     | M A I N   L O O P (Reverse communication loop) |
     %------------------------------------------------% */
  
  while (1) {
    
    /* %---------------------------------------------%
       | Repeatedly call the routine DSAUPD and take | 
       | actions indicated by parameter IDO until    |
       | either convergence is indicated or maxitr   |
       | has been exceeded.                          |
       %---------------------------------------------% */
    
    IGRAPH_F77(dsaupd, DSAUPD) (&ido, bmat, &n, which, &nev, &tol, resid,
				&ncv, v, &ldv, iparam, ipntr, workd, workl,
				&lworkl, &info);
    
    if (ido==-1 || ido==1) {
      
      /* %--------------------------------------%
	 | Perform matrix vector multiplication |
	 |              y <--- OP*x             |
	 | The user should supply his/her own   |
	 | matrix vector multiplication routine |
	 | here that takes workd(ipntr(1)) as   |
	 | the input, and return the result to  |
	 | workd(ipntr(2)).                     |
	 %--------------------------------------% */
      
      int i, j, nlen;
      igraph_real_t *from=&workd[ipntr[0]-1];
      igraph_real_t *to=&workd[ipntr[1]-1];

      for (i=0; i<no_of_nodes; i++) {
	IGRAPH_CHECK(igraph_neighbors(graph, &neis, i, IGRAPH_ALL));
	nlen=igraph_vector_size(&neis);
	to[i]=0.0;
	for (j=0; j<nlen; j++) {
	  long int nei=VECTOR(neis)[j];
	  to[i] += from[nei];
	}
      }				      
      
    } else {
      break;
    }
    
  }

  iparam[4]=1;			/* get eigenvector */
  IGRAPH_F77(dseupd, DSEUPD) (&rvec, all, select, d, v, &ldv, &sigma,
			      bmat, &n, which, &nev, &tol, resid, &ncv,
			      v, &ldv, iparam, ipntr, workd, workl, &lworkl,
			      &ierr);
  
  printf("Eigenvalue: %f\n", d[0]);
  printf("Eigenvector: ");
  for (i=0; i<n; i++) {
    printf("%f ", v[i]);
  }
  printf("\n");
  
  /* Free resources */
  igraph_free(select);
  igraph_free(ax);
  igraph_free(resid);
  igraph_free(d);
  igraph_free(workd);
  igraph_free(workl);
  igraph_free(v);

  igraph_vector_destroy(&neis);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}
