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

SEXP REST_measure_dynamics_idage(SEXP interface, SEXP graph, SEXP st,
				 SEXP pagebins, SEXP pmaxind, SEXP psd) {
  
  REST_i_ptrtable_t ptrtable = REST_i_getptrtable(graph);
  
  SEXP result;
  SEXP akl, sd;
  long int agebins;
  long int maxind;
  long int no_of_nodes;
  long int binwidth;
  
  int *indegree;
  SEXP ntkl, ch, normfact, notnull;
  SEXP tmp, outmode;
  SEXP neis;
  
  long int node;
  long int i, j, k;
  long int edges=0;
  int lsd;

/*   int mptr=1; int m_x=3, m_y=10; */
  
  agebins=R(pagebins);
  maxind=R(pmaxind);
  no_of_nodes=R(ptrtable.vcount(interface, graph));
  indegree=(int*) R_alloc(no_of_nodes, sizeof(int));
  memset(indegree, 0, sizeof(int)*no_of_nodes);
  binwidth=(long int) ceil(((double)(no_of_nodes-1))/agebins);
  lsd=LOGICAL(psd)[0];

  PROTECT(result=NEW_LIST(3));
  SET_VECTOR_ELT(result, 0, NEW_NUMERIC(agebins*(maxind+1)));
  akl=VECTOR_ELT(result, 0);
  if (lsd) { 
    SET_VECTOR_ELT(result, 1, NEW_NUMERIC(agebins*(maxind+1)));
    sd=VECTOR_ELT(result, 1);
  } else {
    SET_VECTOR_ELT(result, 1, NULL_USER_OBJECT);
  }
/*   SET_VECTOR_ELT(result, 2, NEW_NUMERIC(10000)); */
/*   memset(REAL(VECTOR_ELT(result,2)), 0, sizeof(double)*10000); */
  PROTECT(ntkl=NEW_NUMERIC(agebins*(maxind+1)));
  PROTECT(ch=NEW_NUMERIC(agebins*(maxind+1)));
  PROTECT(normfact=NEW_NUMERIC(agebins*(maxind+1)));
  PROTECT(notnull=NEW_NUMERIC(agebins*(maxind+1)));
  memset(REAL(akl), 0, sizeof(double)*agebins*(maxind+1));
  memset(REAL(ntkl), 0, sizeof(double)*agebins*(maxind+1));
  memset(REAL(ch), 0, sizeof(double)*agebins*(maxind+1));
  memset(REAL(normfact), 0, sizeof(double)*agebins*(maxind+1));
  memset(REAL(notnull), 0, sizeof(double)*agebins*(maxind+1));

  PROTECT(tmp=NEW_INTEGER(2));
  INTEGER(tmp)[0]=maxind+1;
  INTEGER(tmp)[1]=agebins;
  SET_DIM(akl, tmp);
  SET_DIM(ntkl, tmp);
  SET_DIM(ch, tmp);
  SET_DIM(normfact, tmp);
  SET_DIM(notnull, tmp);
  if (lsd) {
    SET_DIM(sd, tmp);
  }
  
  PROTECT(outmode=ScalarString(CREATE_STRING_VECTOR("out")));

  for (node=1; node<=no_of_nodes; node++) {
    
    /* inspect the edges */
    neis=ptrtable.neighbors(interface, graph, node, outmode);
    for (i=0; i<GET_LENGTH(neis); i++) {
      long int to=REAL(neis)[i];
      long int xidx=indegree[to-1]+1;
      long int yidx=(long int) ceil(((double)(node-to))/binwidth);
      
      double xk=REAL(st)[node-1]/RMATRIX(ntkl, xidx, yidx);
      double oldm=RMATRIX(akl, xidx, yidx);
      RMATRIX(notnull, xidx, yidx) += 1;
      RMATRIX(akl, xidx, yidx) += (xk-oldm)/RMATRIX(notnull, xidx, yidx);
      if (lsd) {
	RMATRIX(sd, xidx, yidx) += (xk-oldm)*(xk-RMATRIX(akl, xidx, yidx));
      }

/*       if (xidx==m_x && yidx==m_y) { */
/* 	REAL(VECTOR_ELT(result, 2))[mptr++] = xk; */
/*       } */

      indegree[to-1] ++;
      RMATRIX(ntkl, xidx, yidx)--;
      if (RMATRIX(ntkl,xidx, yidx)==0) {
	RMATRIX(normfact, xidx, yidx) += (edges-RMATRIX(ch, xidx, yidx)+1);
	RMATRIX(ch, xidx, yidx)=edges;
      }
      RMATRIX(ntkl, xidx+1, yidx)++;
      if (RMATRIX(ntkl, xidx+1, yidx)==1) {
	RMATRIX(ch, xidx+1, yidx)=edges;
      }
      edges++;
    }

    /* new node, aging */
    RMATRIX(ntkl, 1, 1)++;
    if (RMATRIX(ntkl, 1, 1)==1) {
      RMATRIX(ch, 1, 1)=edges;
    }
    for (k=1; node-1-binwidth*k+1 >=1; k++) {
      long int shnode=node-1-binwidth*k+1;
      long int deg=indegree[shnode-1]+1;
      RMATRIX(ntkl, deg, k)--;
      if (RMATRIX(ntkl, deg, k)==0) {
	RMATRIX(normfact, deg, k) += (edges-RMATRIX(ch, deg, k)+1);
	RMATRIX(ch, deg, k)=edges;
      }
      RMATRIX(ntkl, deg, k+1)++;
      if (RMATRIX(ntkl, deg, k+1)==1) {
	RMATRIX(ch, deg, k+1)=edges;
      }
    }
  }

  /* Ok, measurement done, update change */
  for (i=1; i<=maxind+1; i++) {
    for (j=1; j<=agebins; j++) {
      if (RMATRIX(ntkl, i, j) != 0) {
	RMATRIX(normfact, i, j) += (edges-RMATRIX(ch, i, j)+1);
      }
      RMATRIX(akl, i, j) *= RMATRIX(notnull, i, j) / RMATRIX(normfact, i, j);
      if (lsd) {
	RMATRIX(sd, i, j) += 
	  RMATRIX(akl, i, j)*RMATRIX(akl,i,j)*
	  (RMATRIX(normfact,i,j)/RMATRIX(notnull,i,j)-1);
	if (RMATRIX(normfact,i,j) > 0) {
	  RMATRIX(sd, i, j) = 
	    sqrt(RMATRIX(sd, i, j)/(RMATRIX(normfact,i,j)-1));
	  RMATRIX(sd, i, j) = 
	    2 * RMATRIX(sd,i,j)/sqrt(RMATRIX(normfact,i,j));
	}
      }
    }
  }

/*   REAL(VECTOR_ELT(result,2))[0]=RMATRIX(normfact, m_x, m_y); */
/*   for (i=1; i<mptr; i++) { */
/*     REAL(VECTOR_ELT(result, 2))[i] /= RMATRIX(normfact, m_x, m_y); */
/*   } */
  
  UNPROTECT(7);
  return result;
}

SEXP REST_measure_dynamics_idage_st(SEXP interface, SEXP graph, SEXP akl,
				    SEXP pagebins, SEXP pmaxind) {

  REST_i_ptrtable_t ptrtable = REST_i_getptrtable(graph);

  SEXP result;
  long int agebins;
  long int maxind;
  long int no_of_nodes;
  long int binwidth;
  
  int *indegree;
  SEXP outmode, neis;
  
  long int node;
  long int i, j, k;
  
  agebins=R(pagebins);
  maxind=R(pmaxind);
  no_of_nodes=R(ptrtable.vcount(interface, graph));
  indegree=(int*) R_alloc(no_of_nodes, sizeof(int));
  memset(indegree, 0, sizeof(int)*no_of_nodes);
  binwidth=(long int) ceil(((double)(no_of_nodes-1))/agebins);
  
  PROTECT(result=NEW_NUMERIC(no_of_nodes));
  memset(REAL(result), 0, sizeof(double)*no_of_nodes);
  PROTECT(outmode=ScalarString(CREATE_STRING_VECTOR("out")));
  REAL(result)[0]=RMATRIX(akl, 1, 1);

  for (node=2; node<=no_of_nodes; node++) {
    
    /* new node, aging */
    REAL(result)[node-1] = REAL(result)[node-2] + RMATRIX(akl, 1, 1);
    for (k=1; node-1-binwidth*k+1 >= 1; k++) {
      long int shnode=node-1-binwidth*k+1;
      long int deg=indegree[shnode-1]+1;
      REAL(result)[node-1] += -RMATRIX(akl, deg, k)+RMATRIX(akl, deg, k+1);
    }
    
    /* inspect the outgoing edges */
    neis=ptrtable.neighbors(interface, graph, node, outmode);
    for (i=0; i<GET_LENGTH(neis); i++) {
      long int to=REAL(neis)[i];
      long int xidx=indegree[to-1]+1;
      long int yidx=(long int) ceil(((double)(node-to))/binwidth);
      
      indegree[to-1] ++;
      
      REAL(result)[node-1] += 
	-RMATRIX(akl, xidx, yidx) + RMATRIX(akl, xidx+1, yidx);
    }
  }
  
  UNPROTECT(2);
  return result;
}

