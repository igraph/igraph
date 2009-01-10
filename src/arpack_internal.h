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

#ifndef ARPACK_INTERNAL_H
#define ARPACK_INTERNAL_H

/* Note: only files calling the arpack routines directly need to
   include this header.
*/

#include "igraph.h"
#include "config.h"

#ifndef INTERNAL_ARPACK
#define igrapharpackdsaupd_   dsaupd_
#define igrapharpackdseupd_   dseupd_
#define igrapharpackdnaupd_   dnaupd_
#define igrapharpackdneupd_   dneupd_
#endif

int igrapharpackdsaupd_(long int *ido, char *bmat, long int *n,
		  char *which, long int *nev, igraph_real_t *tol,
		  igraph_real_t *resid, long int *ncv, igraph_real_t *v,
		  long int *ldv, long int *iparam, long int *ipntr, 
		  igraph_real_t *workd, igraph_real_t *workl,
		  long int *lworkl, long int *info);

int igrapharpackdseupd_(long int *rvec, char *howmny, long int *select,
		  igraph_real_t *d, igraph_real_t *z, long int *ldz,
		  igraph_real_t *sigma, char *bmat, long int *n,
		  char *which, long int *nev, igraph_real_t *tol,
		  igraph_real_t *resid, long int *ncv, igraph_real_t *v,
		  long int *ldv, long int *iparam, long int *ipntr, 
		  igraph_real_t *workd, igraph_real_t *workl,
		  long int *lworkl, long int *info);

int igrapharpackdnaupd_(long int *ido, char *bmat, long int *n,
		  char *which, long int *nev, igraph_real_t *tol,
		  igraph_real_t *resid, long int *ncv, igraph_real_t *v,
		  long int *ldv, long int *iparam, long int *ipntr, 
		  igraph_real_t *workd, igraph_real_t *workl,
		  long int *lworkl, long int *info);

int igrapharpackdneupd_(long int *rvec, char *howmny, long int *select,
		  igraph_real_t *dr, igraph_real_t *di,
		  igraph_real_t *z, long int *ldz,
		  igraph_real_t *sigmar, igraph_real_t *sigmai, 
		  igraph_real_t *workev, char *bmat, long int *n,
		  char *which, long int *nev, igraph_real_t *tol,
		  igraph_real_t *resid, long int *ncv, igraph_real_t *v,
		  long int *ldv, long int *iparam, long int *ipntr, 
		  igraph_real_t *workd, igraph_real_t *workl,
		  long int *lworkl, long int *info);
  
#endif
