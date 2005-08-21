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

#define EDGE(i) (REAL(el)[ (long int) REAL(index)[(i)-1] -1 ])

SEXP REST_indexededgelist_create_startindex(SEXP el, SEXP index, SEXP nodes) {
    
  SEXP result, SEXP ;
  long int no_of_nodes;
  long int no_of_edges;
  long int i, j, idx;
  
  no_of_nodes=R(nodes);
  no_of_edges=GET_LENGTH(el);
  
  /* result */
  
  PROTECT(result=NEW_NUMERIC(no_of_nodes+1));
  
  /* create the index */

  idx=0;
  for (i=1; i<=EDGE(1); i++) {
     idx++; REAL(result)[idx-1]=1;
  }
  for (i=1; i<=no_of_edges; i++) {
    long int n=EDGE(i) - EDGE((long int)REAL(result)[idx-1]);
    for (j=0; j<n; j++) {
      idx++; REAL(result)[idx-1]=i;
    }
  }
  j=EDGE((long int)REAL(result)[idx-1]);
  for (i=0; i<=no_of_nodes-j; i++) {
    idx++; REAL(result)[idx-1]=no_of_edges+1;
  }

  /* clean */

  UNPROTECT(1);
  return result;
}
