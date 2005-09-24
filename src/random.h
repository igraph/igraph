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
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/

#ifndef REST_RANDOM_H
#define REST_RANDOM_H

#include "types.h"

#include <stdlib.h>
#include <time.h>

#ifdef USING_R

void GetRNGstate(void);
void PutRNGstate(void);
double  unif_rand(void);
double  norm_rand(void);
double  Rf_rgeom(double);
#define rgeom Rf_rgeom

#define RNG_BEGIN()       GetRNGstate()
#define RNG_END()         PutRNGstate()
#define RNG_INTEGER(l, h) ((long int)(unif_rand()*((h)-(l)+1)+(l)))
#define RNG_NORMAL(m, s)  (norm_rand()*(s)+(m))
#define RNG_UNIF(l, h)    (unif_rand()*((h)-(l))+(l))
#define RNG_UNIF01()      (unif_rand())
#define RNG_GEOM(p)       (rgeom(p))

#else

double igraph_norm_rand(void);
double igraph_rgeom(double);

#define RNG_BEGIN()       srand(time(0))
#define RNG_END()  
#define RNG_INTEGER(l, h) ((long int)((rand())/((double)RAND_MAX+1)*((h)-(l)+1)+(l)))
#define RNG_NORMAL(m, s)  (igraph_norm_rand()*(s)+(m))
#define RNG_UNIF(l, h)    (rand()/((double)RAND_MAX+1)*(double)((h)-(l))+l)
#define RNG_UNIF01()      (RNG_UNIF(0,1))
#define RNG_GEOM(p)       (igraph_rgeom(p))

#endif

#endif
