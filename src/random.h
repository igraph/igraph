/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2003, 2004, 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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

#ifndef REST_RANDOM_H
#define REST_RANDOM_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

#include "types.h"

#include <stdlib.h>
#include <time.h>

#ifdef USING_R

void GetRNGstate(void);
void PutRNGstate(void);
double  unif_rand(void);
double  norm_rand(void);
double  Rf_rgeom(double);
double  Rf_rbinom(double, double);
#define rgeom Rf_rgeom
#define rbinom Rf_rbinom

#define RNG_BEGIN()       GetRNGstate()
#define RNG_END()         PutRNGstate()
#define RNG_INTEGER(l, h) ((long int)(unif_rand()*((h)-(l)+1)+(l)))
#define RNG_NORMAL(m, s)  (norm_rand()*(s)+(m))
#define RNG_UNIF(l, h)    (unif_rand()*((h)-(l))+(l))
#define RNG_UNIF01()      (unif_rand())
#define RNG_GEOM(p)       (rgeom(p))
#define RNG_BINOM(n,p)    (rbinom((n),(p)))
#define RNG_INT31()       ((long)(unif_rand()*0x7FFFFFFFL))

#else

double igraph_norm_rand(void);
double igraph_rgeom(double);
double igraph_rbinom(double, double);
extern int igraph_rng_inited;

#define RNG_BEGIN()       if (!igraph_rng_inited) { srand(time(0)); igraph_rng_inited=1; }
#define RNG_END()  
#define RNG_INTEGER(l, h) ((long int)((rand())/((double)RAND_MAX+1)*((h)-(l)+1)+(l)))
#define RNG_NORMAL(m, s)  (igraph_norm_rand()*(s)+(m))
#define RNG_UNIF(l, h)    (rand()/((double)RAND_MAX+1)*(double)((h)-(l))+(l))
#define RNG_UNIF01()      (RNG_UNIF(0,1))
#define RNG_GEOM(p)       (igraph_rgeom(p))
#define RNG_BINOM(n,p)    (igraph_rbinom((n),(p)))
/* #define RNG_INT31()       ((long)(RNG_UNIF01()*0x7FFFFFFF)) */
#define RNG_INT31()       (rand())

#endif

__END_DECLS

#endif
