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

SEXP REST_cocitation(SEXP interface, SEXP graph, SEXP mode) {

  REST_i_ptrtable_t ptrtable = REST_i_getptrtable(graph);

  SEXP result;
  SEXP dim;

  long int no_of_nodes;
  
  long int from, i, j;
  
  no_of_nodes=R(ptrtable.vcount(interface, graph));
  
  /* The result */
  
  PROTECT(result=NEW_NUMERIC(no_of_nodes*no_of_nodes));
  PROTECT(dim=NEW_INTEGER(2));
  INTEGER(dim)[0]=no_of_nodes;
  INTEGER(dim)[1]=no_of_nodes;
  SET_DIM(result, dim);
  memset(REAL(result), 0, no_of_nodes*no_of_nodes*sizeof(double));
  
  for (from=0; from<no_of_nodes; from++) {
    SEXP neis=ptrtable.neighbors(interface, graph, from+1, mode);
    for (i=0; i < GET_LENGTH(neis)-1; i++) {
      for (j=i+1; j<GET_LENGTH(neis); j++) {
	RMATRIX(result, (long int)REAL(neis)[i], (long int)REAL(neis)[j]) += 1;
	RMATRIX(result, (long int)REAL(neis)[j], (long int)REAL(neis)[i]) += 1;
      }
    }
  }
  
  /* Clean up */
  
  UNPROTECT(2);
  return result;
}
 
  
