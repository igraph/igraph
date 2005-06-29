/* -*- mode: C -*-  */
/* 
   IGraph R package.
   Copyright (C) 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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

SEXP REST_running_mean(SEXP data, SEXP pbinwidth) {

  long int binwidth;
  double sum=0;
  SEXP result;
  long int i;

  binwidth=R(pbinwidth);
  
  /* Memory for result */ 

  PROTECT(result=NEW_NUMERIC(GET_LENGTH(data)-binwidth+1));
  
  /* Initial bin */
  for (i=0; i<binwidth; i++) {
    sum += REAL(data)[i];
  }
  
  REAL(result)[0]=sum/binwidth;
  
  for (i=1; i<GET_LENGTH(data)-binwidth+1; i++) {
    sum -= REAL(data)[i];
    sum += REAL(data)[i+binwidth-1];
    REAL(result)[i] = sum/binwidth;
  }
  
  UNPROTECT(1);
  return result;
}
