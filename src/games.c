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

SEXP REST_ba_game(SEXP pn, SEXP pm, SEXP outseq, SEXP poutpref) {
  
  SEXP result;
  
  long int no_of_nodes;
  long int no_of_neighbors;
  int outpref;
  long int *bag;
  long int bagp;
  long int no_of_edges;
  
  long int resp=0;

  long int i,j;

  no_of_nodes=R(pn);
  if (!isNull(pm)) { no_of_neighbors=R(pm); }
  outpref=LOGICAL(poutpref)[0];
  
  if (isNull(outseq)) {
    bag=(long int*) R_alloc(no_of_nodes * no_of_neighbors + no_of_nodes +
			    outpref * no_of_nodes * no_of_neighbors, 
			    sizeof(long int));
    no_of_edges=(no_of_nodes-1)*no_of_neighbors;
  } else {
    no_of_edges=0;
    for (i=1; i<GET_LENGTH(outseq); i++) {
       no_of_edges+=REAL(outseq)[i];
    }
    bag=(long int*) R_alloc(no_of_nodes + no_of_edges + outpref * no_of_edges,
			    sizeof(long int));
  }
  
  PROTECT(result=NEW_NUMERIC(no_of_edges*2));
  
  /* The first node */

  bagp=0;
  bag[bagp++]=0;
  
  GetRNGstate();

  /* and the others */
  
  for (i=1; i<no_of_nodes; i++) {
    /* draw edges */
    if (isNull(pm)) { no_of_neighbors=(long int) (REAL(outseq)[i]); }
    for (j=0; j<no_of_neighbors; j++) {
      long int to=bag[RNG_INTEGER(0, bagp-1)];
      REAL(result)[resp++] = i+1;
      REAL(result)[resp++] = to+1;
    }
    /* update bag */
    bag[bagp++] = i;
    for (j=0; j<no_of_neighbors; j++) {
      bag[bagp++] = (long int) REAL(result)[resp-2*j-1]-1;
      if (outpref) {
	bag[bagp++] = i;
      }
    }
    R_CheckUserInterrupt();
  }

  PutRNGstate();

  /* This is it, clean */
  UNPROTECT(1);
  return result;
}

/**
 */

