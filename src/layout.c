/* -*- mode: C -*-  */
/* 
   IGraph R package.
   Copyright (C) 2003, 2004  Gabor Csardi <csardi@rmki.kfki.hu>
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
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/

#include "igraph.h"

int igraph_layout_random(igraph_t *graph, vector_t *res) {
  /* TODO */
  return 0;
}

int igraph_layout_circle(igraph_t *graph, vector_t *res) {
  /* TODO */
  return 0;
}

int igraph_layout_fruchterman_reingold(igraph_t *graph, vector_t *res, 
				       integer_t niter, real_t coolexp,
				       integer_t frame, vector_t *initial,
				       real_t initemp) {
  /* TODO */
  return 0;
}

int igraph_layout_kamada_kawai(igraph_t *graph, vector_t *res,
			       integer_t niter, real_t sigma, 
			       real_t initemp, real_t coolexp,
			       real_t kkconst) {
  /* TODO */
  return 0;
}

int igraph_layout_springs(igraph_t *graph, vector_t *res,
			  real_t mass, real_t equil, real_t k,
			  real_t repeqdis, real_t kfr, bool_t repulse) {
  /* TODO */
  return 0;
}

/* #include <Rmath.h> */

/* SEXP REST_layout_kamadakawai(SEXP pn, SEXP pniter,  */
/* 			     SEXP pelen, SEXP pinitemp, SEXP pcoolexp,  */
/* 			     SEXP pkkconst, SEXP psigma, SEXP px, SEXP py) { */

/*   SEXP result, dim; */

/*   double initemp, coolexp, sigma, temp, candx, candy; */
/*   double dpot, odis, ndis, osqd, nsqd, kkconst; */
/*   int niter; */
/*   long int n,i,j,k; */
/*   double *x, *y, *elen; */

/*   /\*Define various things*\/ */
/*   n=R(pn); */
/*   niter=I(pniter); */
/*   initemp=R(pinitemp); */
/*   coolexp=R(pcoolexp); */
/*   kkconst=R(pkkconst); */
/*   sigma=R(psigma); */
/*   x=REAL(px); */
/*   y=REAL(py); */
/*   elen=REAL(pelen); */
/*   GetRNGstate();   /\*Get the RNG state*\/ */

/*   PROTECT(result=NEW_NUMERIC(n*2)); */
/*   PROTECT(dim=NEW_INTEGER(2)); */
/*   INTEGER(dim)[0]=n; */
/*   INTEGER(dim)[1]=2; */
/*   SET_DIM(result, dim); */

/*   /\*Perform the annealing loop*\/ */
/*   temp=initemp; */
/*   for(i=0;i<niter;i++){ */
/*     /\*Update each vertex*\/ */
/*     for(j=0;j<n;j++){ */
/*       /\*Draw the candidate via a gaussian perturbation*\/ */
/*       candx=rnorm(x[j],sigma*temp/initemp); */
/*       candy=rnorm(y[j],sigma*temp/initemp); */
/*       /\*Calculate the potential difference for the new position*\/ */
/*       dpot=0.0; */
/*       for(k=0;k<n;k++)  /\*Potential differences for pairwise effects*\/ */
/*         if(j!=k){ */
/*           odis=sqrt((x[j]-x[k])*(x[j]-x[k])+(y[j]-y[k])*(y[j]-y[k])); */
/*           ndis=sqrt((candx-x[k])*(candx-x[k])+(candy-y[k])*(candy-y[k])); */
/*           osqd=(odis-elen[j+k*n])*(odis-elen[j+k*n]); */
/*           nsqd=(ndis-elen[j+k*n])*(ndis-elen[j+k*n]); */
/*           dpot+=kkconst*(osqd-nsqd)/(elen[j+k*n]*elen[j+k*n]); */
/*         } */
/*       /\*Make a keep/reject decision*\/ */
/*       if(log(runif(0.0,1.0))<dpot/temp){ */
/*         REAL(result)[j  ]=x[j]=candx; */
/*         REAL(result)[j+n]=y[j]=candy; */
/*       } */
/*     } */
/*     /\*Cool the system*\/ */
/*     temp*=coolexp; */
/*   } */
/*   PutRNGstate();   /\*Update the RNG*\/ */

/*   UNPROTECT(2); */
/*   return result; */
/* } */
