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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include "igraph.h"
#include "memory.h"
#include "config.h"

#include <math.h>

int igraph_measure_dynamics_id(const igraph_t *graph,
			       igraph_matrix_t *ak, igraph_matrix_t *sd,
			       igraph_matrix_t *no, igraph_vector_t *cites,
			       igraph_vector_t *debug, 
			       igraph_integer_t debugdeg,
			       const igraph_vector_t *st, 
			       igraph_integer_t pmaxind) {
  
  long int maxind=pmaxind;
  long int no_of_nodes=igraph_vcount(graph);

  int *indegree;
  igraph_matrix_t normfact;
  igraph_vector_t ntk, ch, notnull;
  igraph_vector_t neis;
  
  long int node;
  long int i;
  long int edges=0;
  
  igraph_bool_t lsd=(sd != 0);

  igraph_vector_init(&neis, 0);
  indegree=igraph_Calloc(no_of_nodes, int);
  igraph_matrix_resize(ak, maxind+1, 1);
  igraph_matrix_null(ak);
  if (lsd) {
    igraph_matrix_resize(sd, maxind+1, 1);
    igraph_matrix_null(sd);
  }
  igraph_vector_init(&ntk, maxind+1);
  igraph_vector_init(&ch, maxind+1);
  igraph_matrix_init(&normfact, maxind+1, 1);
  igraph_vector_init(&notnull, maxind+1);

  for (node=0; node<no_of_nodes; node++) {

    IGRAPH_ALLOW_INTERRUPTION();
    
    /* estimate Ak */    
    igraph_neighbors(graph, &neis, node, IGRAPH_OUT);
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=indegree[to];
      
      double xk=VECTOR(*st)[node]/VECTOR(ntk)[xidx];
      double oldm=MATRIX(*ak, xidx, 0);
      VECTOR(notnull)[xidx] += 1;
      MATRIX(*ak, xidx, 0) += (xk-oldm)/VECTOR(notnull)[xidx];
      if (lsd) {
	MATRIX(*sd, xidx, 0) += (xk-oldm)*(xk-MATRIX(*ak, xidx, 0));
      }
      if (debug && xidx==debugdeg) {
	igraph_vector_push_back(debug, xk);
      }
    }

    /* Update ntk, ch, normfact */
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=indegree[to];
      
      indegree[to]++;
      VECTOR(ntk)[xidx] --;
      if (VECTOR(ntk)[xidx]==0) {
	MATRIX(normfact, xidx, 0) += (edges-VECTOR(ch)[xidx]);
	VECTOR(ch)[xidx]=edges;
      }
      VECTOR(ntk)[xidx+1]++;
      if (VECTOR(ntk)[xidx+1]==1) {
	VECTOR(ch)[xidx+1]=edges;
      }
      edges++;
    }
    
    VECTOR(ntk)[0]++;
    if (VECTOR(ntk)[0]==1) {
      VECTOR(ch)[0]=edges;
    }
  }      

  /* Ok, measurement done, update change */
  for (i=0; i<maxind+1; i++) {
    igraph_real_t oldmean;
    if (VECTOR(ntk)[i] != 0) {
      MATRIX(normfact, i, 0) += (edges-VECTOR(ch)[i]);
    }
    oldmean=MATRIX(*ak, i, 0);
    MATRIX(*ak, i, 0) *= VECTOR(notnull)[i] / MATRIX(normfact, i, 0);
    if (lsd) {
      MATRIX(*sd, i, 0) +=
	oldmean * oldmean * VECTOR(notnull)[i] *
	(1-VECTOR(notnull)[i]/MATRIX(normfact, i, 0));
      if (MATRIX(normfact, i, 0) > 0) {
	MATRIX(*sd, i, 0) =
	  sqrt(MATRIX(*sd, i, 0)/(MATRIX(normfact, i, 0)-1));
      }
    }
  }

  if (no) {
    igraph_matrix_destroy(no);
    *no=normfact;
  } else {
    igraph_matrix_destroy(&normfact);
  }
  if (cites) {
    igraph_vector_destroy(cites);
    *cites=notnull;
  } else {
    igraph_vector_destroy(&notnull);
  }
  
  igraph_Free(indegree);
  igraph_vector_destroy(&ntk);
  igraph_vector_destroy(&ch);
  igraph_vector_destroy(&neis);
  
  return 0;
}

int igraph_measure_dynamics_id_expected(const igraph_t *graph,
					igraph_vector_t *res,
					const igraph_vector_t *ak,
					const igraph_vector_t *st,
					igraph_integer_t pmaxind) {
  long int maxind=pmaxind;
  
  igraph_vector_t ntk;
  igraph_vector_t indegree;
  igraph_vector_t neis;
  
  long int no_of_nodes=igraph_vcount(graph);
  long int node, i, j;

  IGRAPH_VECTOR_INIT_FINALLY(&ntk, maxind+1);
  IGRAPH_VECTOR_INIT_FINALLY(&indegree, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  
  IGRAPH_CHECK(igraph_vector_resize(res, maxind+1));
  igraph_vector_null(res);

  for (node=0; node<no_of_nodes; node++) {    
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    /* measure expected number of citations */
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, node, IGRAPH_OUT));
    for (i=0; i<igraph_vector_size(&neis); i++) {
      for (j=0; j<maxind+1; j++) {
	VECTOR(*res)[j] +=
	  VECTOR(*ak)[j]*VECTOR(ntk)[j]/VECTOR(*st)[node]; 
      }      
    }

    /* update degree & ntk */
    VECTOR(ntk)[0]++;
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=VECTOR(indegree)[to];
      VECTOR(indegree)[to]++;
      VECTOR(ntk)[xidx]--;
      VECTOR(ntk)[xidx+1]++;
    }
  }

  igraph_vector_destroy(&neis);
  igraph_vector_destroy(&indegree);
  igraph_vector_destroy(&ntk);
  IGRAPH_FINALLY_CLEAN(3);
  
  return 0;
}

int igraph_measure_dynamics_id_expected2(const igraph_t *graph,
					 igraph_vector_t *res,
					 const igraph_vector_t *ak,
					 const igraph_vector_t *st,
					 igraph_integer_t pmaxind) {
  long int maxind=pmaxind;
  
  igraph_vector_t ntk;
  igraph_vector_t cumst;
  igraph_vector_t ch;
  igraph_vector_t indegree;
  igraph_vector_t outdegree;
  igraph_vector_t neis;

  long int no_of_nodes=igraph_vcount(graph);
  long int node, i;
  
  IGRAPH_VECTOR_INIT_FINALLY(&ntk, maxind+1);
  IGRAPH_VECTOR_INIT_FINALLY(&ch, maxind+1);
  IGRAPH_VECTOR_INIT_FINALLY(&cumst, no_of_nodes+1);
  IGRAPH_VECTOR_INIT_FINALLY(&indegree, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&outdegree, no_of_nodes);

  IGRAPH_CHECK(igraph_degree(graph, &outdegree, igraph_vss_all(),
			     IGRAPH_OUT, IGRAPH_LOOPS));

  /* create the cumulative sum of dt/S(t) */
  VECTOR(cumst)[0]=0;
  for (i=0; i<no_of_nodes; i++) {
    VECTOR(cumst)[i+1] = VECTOR(cumst)[i] + 
      VECTOR(outdegree)[i]/VECTOR(*st)[i];
  }

  igraph_vector_destroy(&outdegree);
  IGRAPH_FINALLY_CLEAN(1);

  IGRAPH_CHECK(igraph_vector_resize(res, maxind+1));
  igraph_vector_null(res);
  
  for (node=0; node<no_of_nodes; node++) {
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, node, IGRAPH_OUT));
    
    /* update degree and ntk */
    /* update result if needed */
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=VECTOR(indegree)[to];
      VECTOR(indegree)[to]++;
      
      VECTOR(ntk)[xidx]--;
      VECTOR(*res)[xidx] += (VECTOR(ntk)[xidx]+1)*
	(VECTOR(cumst)[node]-VECTOR(cumst)[(long int)VECTOR(ch)[xidx]]);
      VECTOR(ch)[xidx]=node;

      VECTOR(ntk)[xidx+1]++;
      VECTOR(*res)[xidx+1] += (VECTOR(ntk)[xidx+1]-1)*
	(VECTOR(cumst)[node]-VECTOR(cumst)[(long int)VECTOR(ch)[xidx+1]]);
      VECTOR(ch)[xidx+1]=node;
    }
    
    VECTOR(ntk)[0]++;
    VECTOR(*res)[0] += (VECTOR(ntk)[0]-1)*
      (VECTOR(cumst)[node]-VECTOR(cumst)[(long int)VECTOR(ch)[0]]);
    VECTOR(ch)[0]=node;
  }
  
  /* complete res */
  for (i=0; i<maxind+1; i++) {
    VECTOR(*res)[i] += VECTOR(ntk)[i]*
      (VECTOR(cumst)[node]-VECTOR(cumst)[(long int)VECTOR(ch)[i]]);
    VECTOR(*res)[i] *= VECTOR(*ak)[i];
  }
  
  igraph_vector_destroy(&neis);
  igraph_vector_destroy(&indegree);
  igraph_vector_destroy(&cumst);
  igraph_vector_destroy(&ch);
  igraph_vector_destroy(&ntk);
  IGRAPH_FINALLY_CLEAN(5);
  
  return 0;
}

int igraph_measure_dynamics_id_st(const igraph_t *graph, 
				  igraph_vector_t *res, 
				  const igraph_matrix_t *ak) {

  long int no_of_nodes=igraph_vcount(graph);
  int *indegree;
  igraph_vector_t neis;
  
  long int node;
  long int i;
  
  igraph_vector_init(&neis, 0);
  
  indegree=igraph_Calloc(no_of_nodes, int);
  
  igraph_vector_resize(res, no_of_nodes);
  igraph_vector_null(res);
  VECTOR(*res)[0]=MATRIX(*ak, 0, 0);
  
  for (node=1; node<no_of_nodes; node++) {
    
    /* new node */
    VECTOR(*res)[node] = VECTOR(*res)[node-1] + MATRIX(*ak, 0, 0);
    
    /* outgoing edges */
    igraph_neighbors(graph, &neis, node, IGRAPH_OUT);
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=indegree[to];
      
      indegree[to]++;
      
      VECTOR(*res)[node] += -MATRIX(*ak, xidx, 0)+MATRIX(*ak, xidx+1, 0);
    }
  }
  
  igraph_vector_destroy(&neis);
  igraph_Free(indegree);

  return 0;
}

int igraph_measure_dynamics_idwindow(const igraph_t *graph, 
				     igraph_matrix_t *ak, 
				     igraph_matrix_t *sd,
				     const igraph_vector_t *st,
				     igraph_integer_t pmaxind,
				     igraph_integer_t time_window) {
  long int maxind=pmaxind;
  long int no_of_nodes=igraph_vcount(graph);
  
  igraph_vector_t indegree;
  igraph_matrix_t normfact;
  igraph_dqueue_t history;
  igraph_vector_t ntk, ch, notnull;
  igraph_vector_t neis;
  
  long int node;
  long int i, j;
  long int edges=0;
  
  igraph_bool_t lsd=(sd != 0);
  
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&indegree, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&ntk, maxind+1);
  IGRAPH_VECTOR_INIT_FINALLY(&ch, maxind+1);
  IGRAPH_MATRIX_INIT_FINALLY(&normfact, maxind+1, 1);
  IGRAPH_VECTOR_INIT_FINALLY(&notnull, maxind+1);
  IGRAPH_DQUEUE_INIT_FINALLY(&history, time_window);
  
  igraph_matrix_resize(ak, maxind+1, 1);
  igraph_matrix_null(ak);
  if (lsd) {
    igraph_matrix_resize(sd, maxind+1, 1);
    igraph_matrix_null(sd);
  } 
  
  for (node=0; node<no_of_nodes; node++) {

    IGRAPH_ALLOW_INTERRUPTION();

    /* estimate Ak */
    igraph_neighbors(graph, &neis, node, IGRAPH_OUT);
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=VECTOR(indegree)[to];
      
      double xk=VECTOR(*st)[node]/VECTOR(ntk)[xidx];
      double oldm=MATRIX(*ak, xidx, 0);
      VECTOR(notnull)[xidx]+=1;
      MATRIX(*ak, xidx, 0) += (xk-oldm)/VECTOR(notnull)[xidx];
      if (lsd) {
	MATRIX(*sd, xidx, 0) += (xk-oldm)*(xk-MATRIX(*ak, xidx, 0));
      }
    }

    /* Update ntk, ch, normfact */
    edges += igraph_vector_size(&neis);
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=VECTOR(indegree)[to];
      
      VECTOR(indegree)[to]++;
      VECTOR(ntk)[xidx] --;
      if (VECTOR(ntk)[xidx]==0) {
	MATRIX(normfact, xidx, 0) += (edges-VECTOR(ch)[xidx]);
	VECTOR(ch)[xidx]=edges;
      }
      VECTOR(ntk)[xidx+1]++;
      if (VECTOR(ntk)[xidx+1]==1) {
	VECTOR(ch)[xidx+1]=edges;
      }
      igraph_dqueue_push(&history, to);      
    }
    igraph_dqueue_push(&history, -1);
    
    /* time window */
    if (node > time_window) {
      while ( (j=igraph_dqueue_pop(&history)) != -1) {
	long int xidx=VECTOR(indegree)[j];
	VECTOR(indegree)[j] --; 
	VECTOR(ntk)[xidx]--;
	if (VECTOR(ntk)[xidx]==0) {
	  MATRIX(normfact, xidx, 0) += (edges-VECTOR(ch)[xidx]);
	  VECTOR(ch)[xidx]=edges;
	}
	VECTOR(ntk)[xidx-1]++;
	if (VECTOR(ntk)[xidx-1]==1) {
	  VECTOR(ch)[xidx-1]=edges;
	}
      }
    }

    /* isolate node */
    VECTOR(ntk)[0]++;
    if (VECTOR(ntk)[0]==1) {
      VECTOR(ch)[0]=edges;
    }        
    
  }      

  /* Ok, measurement done, update change */
  for (i=0; i<maxind+1; i++) {
    igraph_real_t oldmean;
    if (VECTOR(ntk)[i] != 0) {
      MATRIX(normfact, i, 0) += (edges-VECTOR(ch)[i]);
    }
    oldmean=MATRIX(*ak, i, 0);
    MATRIX(*ak, i, 0) *= VECTOR(notnull)[i] / MATRIX(normfact, i, 0);
    if (lsd) {
      /* TODO: confidence interval estimation */
      MATRIX(*sd, i, 0) +=
	oldmean * oldmean * VECTOR(notnull)[i] *
	(1-VECTOR(notnull)[i]/MATRIX(normfact, i, 0));
      if (MATRIX(normfact, i, 0) > 0) {
	MATRIX(*sd, i, 0) =
	  sqrt(MATRIX(*sd, i, 0)/(MATRIX(normfact, i, 0)-1));
      }
    }
  }
  
  igraph_dqueue_destroy(&history);
  igraph_vector_destroy(&notnull);
  igraph_matrix_destroy(&normfact);
  igraph_vector_destroy(&ch);
  igraph_vector_destroy(&ntk);
  igraph_vector_destroy(&indegree);
  igraph_vector_destroy(&neis);
  IGRAPH_FINALLY_CLEAN(7);
  
  return 0;
}	   

int igraph_measure_dynamics_idwindow_st(const igraph_t *graph,
					igraph_vector_t *res,
					const igraph_matrix_t *ak,
					igraph_integer_t time_window) {
  long int no_of_nodes=igraph_vcount(graph);
  igraph_vector_t indegree;
  igraph_vector_t neis;
  igraph_dqueue_t history;
  
  long int node;
  long int i, k;
  
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_DQUEUE_INIT_FINALLY(&history, time_window);
  IGRAPH_VECTOR_INIT_FINALLY(&indegree, no_of_nodes);
  
  IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
  igraph_vector_null(res);
  VECTOR(*res)[0]=MATRIX(*ak, 0, 0);
  
  for (node=1; node<no_of_nodes; node++) {
    
    IGRAPH_ALLOW_INTERRUPTION();

    /* new node */
    VECTOR(*res)[node] = VECTOR(*res)[node-1] + MATRIX(*ak, 0, 0);
	
    if(node > time_window) {
      while( (k = igraph_dqueue_pop(&history)) != -1) {
	long int xidx=VECTOR(indegree)[k];
	VECTOR(*res)[node] -= MATRIX(*ak, xidx, 0);
	VECTOR(*res)[node] += MATRIX(*ak, xidx-1, 0);
	VECTOR(indegree)[k]--;
      }
    }
    
    /* outgoing edges */
    igraph_neighbors(graph, &neis, node, IGRAPH_OUT);
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=VECTOR(indegree)[to];
      
      VECTOR(indegree)[to]++;
      igraph_dqueue_push(&history, to);
      
      VECTOR(*res)[node] += -MATRIX(*ak, xidx, 0)+MATRIX(*ak, xidx+1, 0);
    }
    igraph_dqueue_push(&history, -1);
  }
  
  igraph_vector_destroy(&neis);
  igraph_dqueue_destroy(&history);
  igraph_vector_destroy(&indegree);
  IGRAPH_FINALLY_CLEAN(3);

  return 0;
}
    
int igraph_measure_dynamics_idage(const igraph_t *graph,
				  igraph_matrix_t *akl, 
				  igraph_matrix_t *sd,
				  igraph_matrix_t *no,
				  igraph_matrix_t *cites,
				  const igraph_vector_t *st, igraph_integer_t pagebins,
				  igraph_integer_t pmaxind) {

  long int agebins=pagebins;
  long int maxind=pmaxind;
  long int no_of_nodes=igraph_vcount(graph);
  long int binwidth;
  
  int *indegree;
  igraph_matrix_t ntkl, ch, normfact, notnull;
  igraph_vector_t neis;
  
  long int node;
  long int i, j, k;
  long int edges=0;

  igraph_bool_t lsd=(sd != 0);

  binwidth = no_of_nodes/agebins+1;

  igraph_vector_init(&neis, 0);
  indegree=igraph_Calloc(no_of_nodes, int);
  igraph_matrix_resize(akl, maxind+1, agebins);
  igraph_matrix_null(akl);
  if (lsd) {
    igraph_matrix_resize(sd, maxind+1, agebins);
    igraph_matrix_null(sd);
  }
  igraph_matrix_init(&ntkl, maxind+1, agebins+1);
  igraph_matrix_init(&ch, maxind+1, agebins+1);
  igraph_matrix_init(&normfact, maxind+1, agebins);
  igraph_matrix_init(&notnull, maxind+1, agebins);
  
  for (node=0; node<no_of_nodes; node++) {
    
    IGRAPH_ALLOW_INTERRUPTION();

    /* inspect the edges, update A(k,l) */
    igraph_neighbors(graph, &neis, node, IGRAPH_OUT);
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=indegree[to];
      long int yidx=(node-to)/binwidth;
      
      double xk=VECTOR(*st)[node]/MATRIX(ntkl, xidx, yidx);
      double oldm=MATRIX(*akl, xidx, yidx);
      MATRIX(notnull, xidx, yidx) += 1;
      MATRIX(*akl, xidx, yidx) += (xk-oldm)/MATRIX(notnull, xidx, yidx);
      if (lsd) {
	MATRIX(*sd, xidx, yidx) += (xk-oldm)*(xk-MATRIX(*akl, xidx, yidx));
      }
    }

    /* add the new edges */
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=indegree[to];
      long int yidx=(node-to)/binwidth;
      
      indegree[to] ++;
      MATRIX(ntkl, xidx, yidx)--;
      if (MATRIX(ntkl,xidx, yidx)==0) {
	MATRIX(normfact, xidx, yidx) += (edges-MATRIX(ch, xidx, yidx)+1);
	MATRIX(ch, xidx, yidx)=edges;
      }
      MATRIX(ntkl, xidx+1, yidx)++;
      if (MATRIX(ntkl, xidx+1, yidx)==1) {
	MATRIX(ch, xidx+1, yidx)=edges;
      }
      edges++;
    }      

    /* new node, aging */
    MATRIX(ntkl, 0, 0)++;
    if (MATRIX(ntkl, 0, 0)==1) {
      MATRIX(ch, 0, 0)=edges;
    }
    for (k=1; node-binwidth*k+1 >=1; k++) {
      long int shnode=node-binwidth*k;
      long int deg=indegree[shnode];
      MATRIX(ntkl, deg, k-1)--;
      if (MATRIX(ntkl, deg, k-1)==0) {
	MATRIX(normfact, deg, k-1) += (edges-MATRIX(ch, deg, k-1)+1);
	MATRIX(ch, deg, k-1)=edges;
      }
      MATRIX(ntkl, deg, k)++;
      if (MATRIX(ntkl, deg, k)==1) {
	MATRIX(ch, deg, k)=edges;
      }
    }
  }

  /* Ok, measurement done, update change */
  for (i=0; i<maxind+1; i++) {
    for (j=0; j<agebins; j++) {
      igraph_real_t oldmean;
      if (MATRIX(ntkl, i, j) != 0) {
	MATRIX(normfact, i, j) += (edges-MATRIX(ch, i, j)+1);
      }
      oldmean=MATRIX(*akl, i, j);
      MATRIX(*akl, i, j) *= MATRIX(notnull, i, j) / MATRIX(normfact, i, j);
      if (lsd) {
	MATRIX(*sd, i, j) +=
	  oldmean * oldmean * MATRIX(notnull, i, j) * 
	  (1-MATRIX(notnull,i,j)/MATRIX(normfact,i,j));
	if (MATRIX(normfact,i,j) > 0) {
	  MATRIX(*sd, i, j) =
	    sqrt(MATRIX(*sd, i, j)/(MATRIX(normfact,i,j)-1));
	}
      }
    }
  }

  if (no) {
    igraph_matrix_destroy(no);
    *no=normfact;
  } else {
    igraph_matrix_destroy(&normfact);
  }
  if (cites) {
    igraph_matrix_destroy(cites);
    *cites=notnull;
  } else {
    igraph_matrix_destroy(&notnull);
  }

  igraph_Free(indegree);
  igraph_matrix_destroy(&ntkl);
  igraph_matrix_destroy(&ch);
  igraph_vector_destroy(&neis);

  return 0;
}

int igraph_measure_dynamics_idage_expected(const igraph_t *graph,
					   igraph_matrix_t *res,
					   const igraph_matrix_t *akl,
					   const igraph_vector_t *st,
					   igraph_integer_t pmaxind) {
  long int agebins=igraph_matrix_ncol(akl);
  long int no_of_nodes=igraph_vcount(graph);
  long int binwidth=no_of_nodes/agebins+1;
  long int maxind=pmaxind;
  
  igraph_vector_t indegree;
  igraph_vector_t neis;
  igraph_matrix_t ntkl;

  long int node, i, j, k;

  IGRAPH_VECTOR_INIT_FINALLY(&indegree, no_of_nodes);
  IGRAPH_MATRIX_INIT_FINALLY(&ntkl, maxind+1, agebins+1);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  
  IGRAPH_CHECK(igraph_matrix_resize(res, maxind+1, agebins));
  igraph_matrix_null(res);
  
  for (node=0; node<no_of_nodes; node++) {
    long int n;
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    /* expected number of citations */
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, node, IGRAPH_OUT));
    n=igraph_vector_size(&neis);
    for (j=0; j<maxind+1; j++) {
      for (k=0; k<agebins; k++) {
	MATRIX(*res, j, k) += 
	  n * MATRIX(*akl, j, k)*MATRIX(ntkl, j, k)/VECTOR(*st)[node];
      }
    }
    
    /* update degree & ntkl */
    MATRIX(ntkl, 0, 0) += 1;
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=VECTOR(indegree)[to];
      long int yidx=(node-to)/binwidth;
      MATRIX(ntkl, xidx, yidx) -= 1;
      MATRIX(ntkl, xidx+1, yidx) += 1;
      VECTOR(indegree)[to] += 1;
    }
    for (k=1; node-binwidth*k+1 >= 1; k++) {
      long int shnode=node-binwidth*k;
      long int deg=VECTOR(indegree)[shnode];
      MATRIX(ntkl, deg, k-1) -= 1;
      MATRIX(ntkl, deg, k) += 1;
    }
  }
  
  igraph_vector_destroy(&neis);
  igraph_matrix_destroy(&ntkl);
  igraph_vector_destroy(&indegree);
  IGRAPH_FINALLY_CLEAN(3);
  
  return 0;
}

int igraph_measure_dynamics_idage_st(const igraph_t *graph, igraph_vector_t *res,
				     const igraph_matrix_t *akl) {

  long int agebins=igraph_matrix_ncol(akl);
  long int no_of_nodes=igraph_vcount(graph);
  long int binwidth;
  
  int *indegree;
  igraph_vector_t neis;
  
  long int node;
  long int i, k;

  igraph_vector_init(&neis, 0);
  
  indegree=igraph_Calloc(no_of_nodes, int);
  binwidth=no_of_nodes/agebins+1;
  
  igraph_vector_resize(res, no_of_nodes);
  igraph_vector_null(res);
  VECTOR(*res)[0]=MATRIX(*akl, 0, 0);

  for (node=1; node<no_of_nodes; node++) {

    IGRAPH_ALLOW_INTERRUPTION();
    
    /* new node, aging */
    VECTOR(*res)[node] = VECTOR(*res)[node-1] + MATRIX(*akl, 0, 0);
    for (k=1; node-binwidth*k+1 >= 1; k++) {
      long int shnode=node-binwidth*k;
      long int deg=indegree[shnode];
      VECTOR(*res)[node] += -MATRIX(*akl, deg, k-1)+MATRIX(*akl, deg, k);
    }
    
    /* inspect the outgoing edges */
    igraph_neighbors(graph, &neis, node, IGRAPH_OUT);
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=indegree[to];
      long int yidx=(node-to)/binwidth;
      
      indegree[to] ++;
      
      VECTOR(*res)[node] +=
	-MATRIX(*akl, xidx, yidx) + MATRIX(*akl, xidx+1, yidx);
    }
  }
  
  igraph_vector_destroy(&neis);
  igraph_Free(indegree);
  
  return 0;
}

int igraph_measure_dynamics_idwindowage(const igraph_t *graph, 
					igraph_matrix_t *akl, 
					igraph_matrix_t *sd, 
					const igraph_vector_t *st, 
					igraph_integer_t pagebins,
					igraph_integer_t pmaxind, 
					igraph_integer_t time_window) {

  long int agebins=pagebins;
  long int maxind=pmaxind;
  long int no_of_nodes=igraph_vcount(graph);
  long int binwidth;
  
  int *indegree;
  igraph_matrix_t ntkl, ch, normfact, notnull;
  igraph_vector_t neis;
  
  long int node;
  long int i, j, k;
  long int edges=0;

  igraph_bool_t lsd=(sd != 0);

  igraph_dqueue_t history;

  binwidth = no_of_nodes/agebins+1;

  igraph_vector_init(&neis, 0);
  indegree=igraph_Calloc(no_of_nodes, int);
  igraph_matrix_resize(akl, maxind+1, agebins);
  igraph_matrix_null(akl);
  if (lsd) {
    igraph_matrix_resize(sd, maxind+1, agebins);
    igraph_matrix_null(sd);
  }
  igraph_matrix_init(&ntkl, maxind+1, agebins+1);
  igraph_matrix_init(&ch, maxind+1, agebins+1);
  igraph_matrix_init(&normfact, maxind+1, agebins);
  igraph_matrix_init(&notnull, maxind+1, agebins);
  igraph_dqueue_init(&history, time_window);
  
  igraph_dqueue_push(&history, -1);

  for (node=0; node<no_of_nodes; node++) {
    
    IGRAPH_ALLOW_INTERRUPTION();

    /* inspect the edges */
   
    igraph_neighbors(graph, &neis, node, IGRAPH_OUT);
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=indegree[to];
      long int yidx=(node-to)/binwidth;
      
      double xk=VECTOR(*st)[node]/MATRIX(ntkl, xidx, yidx);
      double oldm=MATRIX(*akl, xidx, yidx);
      MATRIX(notnull, xidx, yidx) += 1;
      MATRIX(*akl, xidx, yidx) += (xk-oldm)/MATRIX(notnull, xidx, yidx);
      if (lsd) {
	MATRIX(*sd, xidx, yidx) += (xk-oldm)*(xk-MATRIX(*akl, xidx, yidx));
      }

      indegree[to] ++;
      MATRIX(ntkl, xidx, yidx)--;
      if (MATRIX(ntkl,xidx, yidx)==0) {
	MATRIX(normfact, xidx, yidx) += (edges-MATRIX(ch, xidx, yidx)+1);
	MATRIX(ch, xidx, yidx)=edges;
      }
      MATRIX(ntkl, xidx+1, yidx)++;
      if (MATRIX(ntkl, xidx+1, yidx)==1) {
	MATRIX(ch, xidx+1, yidx)=edges;
      }
      edges++;
      igraph_dqueue_push(&history, to);
    }
    igraph_dqueue_push(&history, -1);

    /* new node, aging */
    MATRIX(ntkl, 0, 0)++;
    if (MATRIX(ntkl, 0, 0)==1) {
      MATRIX(ch, 0, 0)=edges;
    }
    for (k=1; node-binwidth*k+1 >=1; k++) {
      long int shnode=node-binwidth*k;
      long int deg=indegree[shnode];
      MATRIX(ntkl, deg, k-1)--;
      if (MATRIX(ntkl, deg, k-1)==0) {
	MATRIX(normfact, deg, k-1) += (edges-MATRIX(ch, deg, k-1)+1);
	MATRIX(ch, deg, k-1)=edges;
      }
      MATRIX(ntkl, deg, k)++;
      if (MATRIX(ntkl, deg, k)==1) {
	MATRIX(ch, deg, k)=edges;
      }
    }
    
    /* time window */
    if (node > time_window) {
      while ( (j=igraph_dqueue_pop(&history)) != -1) {
	long int xidx=indegree[j];
	long int yidx=(node-j)/binwidth;
	indegree[j]--;
	MATRIX(ntkl, xidx, yidx)--;
	if (MATRIX(ntkl, xidx, yidx)==0) {
	  MATRIX(normfact, xidx, yidx) += (edges-MATRIX(ch, xidx, yidx)+1);
	  MATRIX(ch, xidx, yidx)=edges;
	}
	MATRIX(ntkl, xidx-1, yidx)++;
	if (MATRIX(ntkl, xidx-1, yidx)==1) {
	  MATRIX(ch, xidx-1, yidx)=edges;
	}
      }
    }
  }

  /* Ok, measurement done, update change */
  for (i=0; i<maxind+1; i++) {
    for (j=0; j<agebins; j++) {
      igraph_real_t oldmean;
      if (MATRIX(ntkl, i, j) != 0) {
	MATRIX(normfact, i, j) += (edges-MATRIX(ch, i, j)+1);
      }
      oldmean=MATRIX(*akl, i, j);
      MATRIX(*akl, i, j) *= MATRIX(notnull, i, j) / MATRIX(normfact, i, j);
      if (lsd) {
	MATRIX(*sd, i, j) +=
	  oldmean * oldmean * MATRIX(notnull, i, j) * 
	  (1-MATRIX(notnull,i,j)/MATRIX(normfact,i,j));
	if (MATRIX(normfact,i,j) > 0) {
	  MATRIX(*sd, i, j) =
	    sqrt(MATRIX(*sd, i, j)/(MATRIX(normfact,i,j)-1));
	}
      }
    }
  }
  
  igraph_matrix_destroy(&normfact);

  igraph_dqueue_destroy(&history);
  igraph_Free(indegree);
  igraph_matrix_destroy(&ntkl);
  igraph_matrix_destroy(&ch);
  igraph_matrix_destroy(&notnull);
  igraph_vector_destroy(&neis);

  return 0;
}

int igraph_measure_dynamics_idwindowage_st(const igraph_t *graph, 
					   igraph_vector_t *res,
					   const igraph_matrix_t *akl,
					   igraph_integer_t time_window) {

  long int agebins=igraph_matrix_ncol(akl);
  long int no_of_nodes=igraph_vcount(graph);
  long int binwidth;
  
  int *indegree;
  igraph_vector_t neis;
  
  long int node;
  long int i, k;

  igraph_dqueue_t history;

  igraph_vector_init(&neis, 0);
  igraph_dqueue_init(&history, time_window);
  
  indegree=igraph_Calloc(no_of_nodes, int);
  binwidth=no_of_nodes/agebins+1;
  
  igraph_vector_resize(res, no_of_nodes);
  igraph_vector_null(res);
  VECTOR(*res)[0]=MATRIX(*akl, 0, 0);

  for (node=1; node<no_of_nodes; node++) {

    IGRAPH_ALLOW_INTERRUPTION();
    
    /* new node, aging */
    VECTOR(*res)[node] = VECTOR(*res)[node-1] + MATRIX(*akl, 0, 0);
    for (k=1; node-binwidth*k+1 >= 1; k++) {
      long int shnode=node-binwidth*k;
      long int deg=indegree[shnode];
      VECTOR(*res)[node] += -MATRIX(*akl, deg, k-1)+MATRIX(*akl, deg, k);
    }

    if (node > time_window) {
      while ( (k=igraph_dqueue_pop(&history)) != -1) {
	long int xidx=indegree[k];
	long int yidx=(node-k)/binwidth;
	VECTOR(*res)[node] -= MATRIX(*akl, xidx, yidx);
	VECTOR(*res)[node] += MATRIX(*akl, xidx-1, yidx);
	indegree[k]--;
      }
    }
    
    /* inspect the outgoing edges */
    igraph_neighbors(graph, &neis, node, IGRAPH_OUT);
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=indegree[to];
      long int yidx=(node-to)/binwidth;
      
      indegree[to] ++;
      igraph_dqueue_push(&history, to);
      
      VECTOR(*res)[node] +=
	-MATRIX(*akl, xidx, yidx) + MATRIX(*akl, xidx+1, yidx);
    }
    igraph_dqueue_push(&history, -1);
  }
  
  igraph_vector_destroy(&neis);
  igraph_Free(indegree);
  
  return 0;
}

int igraph_measure_dynamics_citedcat_id_age(const igraph_t *graph,
					    igraph_array3_t *adkl,
					    igraph_array3_t *sd,
					    const igraph_vector_t *st,
					    const igraph_vector_t *cats,
					    igraph_integer_t pno_cats,
					    igraph_integer_t pagebins,
					    igraph_integer_t pmaxind) {
  
  long int agebins=pagebins;
  long int maxind=pmaxind;
  long int no_cats=pno_cats;
  long int no_of_nodes=igraph_vcount(graph);
  long int binwidth;
  
  int *indegree;
  igraph_array3_t ntkl, ch, normfact, notnull;
  igraph_vector_t neis;
  
  long int node;
  long int i,j,k;
  long int edges=0;
  
  igraph_bool_t lsd=(sd != 0);

  binwidth=no_of_nodes/agebins+1;
  
  igraph_vector_init(&neis, 0);
  indegree=igraph_Calloc(no_of_nodes, int);
  igraph_array3_resize(adkl, no_cats, maxind+1, agebins);
  igraph_array3_null(adkl);
  if (lsd) {
    igraph_array3_resize(sd, no_cats, maxind+1, agebins);
    igraph_array3_null(sd);
  }
  igraph_array3_init(&ntkl,     no_cats, maxind+1, agebins);
  igraph_array3_init(&ch,       no_cats, maxind+1, agebins);
  igraph_array3_init(&normfact, no_cats, maxind+1, agebins);
  igraph_array3_init(&notnull,  no_cats, maxind+1, agebins);
  
  for (node=0; node < no_of_nodes; node++) {
    long int cidx;
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    /* update A(d,k,l) */
    igraph_neighbors(graph, &neis, node, IGRAPH_OUT);
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int cidx=VECTOR(*cats)[to];
      long int xidx=indegree[to];
      long int yidx=(node-to)/binwidth;
      
      double xk=VECTOR(*st)[node] / ARRAY3(ntkl, cidx, xidx, yidx);
      double oldm=ARRAY3(*adkl, cidx, xidx, yidx);
      ARRAY3(notnull, cidx, xidx, yidx) += 1;
      ARRAY3(*adkl, cidx, xidx, yidx) += 
	(xk-oldm)/ARRAY3(notnull, cidx, xidx, yidx);
      if (lsd) {
	ARRAY3(*sd, cidx, xidx, yidx) += 
	  (xk-oldm)*(xk-ARRAY3(*adkl, cidx, xidx, yidx));
      }
    }
    
    /* add the new edges */
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int cidx=VECTOR(*cats)[to];
      long int xidx=indegree[to];
      long int yidx=(node-to)/binwidth;
      
      indegree[to] ++;
      ARRAY3(ntkl, cidx, xidx, yidx)--;
      if (ARRAY3(ntkl, cidx, xidx, yidx)==0) {
	ARRAY3(normfact, cidx, xidx, yidx) += 
	  (edges-ARRAY3(ch, cidx, xidx, yidx)+1);
	ARRAY3(ch, cidx, xidx, yidx)=edges;
      }
      ARRAY3(ntkl, cidx, xidx+1, yidx)++;
      if (ARRAY3(ntkl, cidx, xidx+1, yidx)==1) {
	ARRAY3(ch, cidx, xidx+1, yidx)=edges;
      }
      edges++;
    }
    
    /* aging */
    cidx=VECTOR(*cats)[node];
    ARRAY3(ntkl, cidx, 0, 0)++;
    if (ARRAY3(ntkl, cidx, 0, 0)==1) {
      ARRAY3(ch, cidx, 0, 0)=edges;
    }
    for (k=1; node-binwidth*k+1 >= 1; k++) {
      long int shnode=node-binwidth*k;
      long int cat=VECTOR(*cats)[shnode];
      long int deg=indegree[shnode];
      ARRAY3(ntkl, cat, deg, k-1)--;
      if (ARRAY3(ntkl, cat, deg, k-1)==0) {
	ARRAY3(normfact, cat, deg, k-1)+=(edges-ARRAY3(ch, cat, deg, k-1)+1);
	ARRAY3(ch, cat, deg, k-1)=edges;
      }
      ARRAY3(ntkl, cat, deg, k)++;
      if (ARRAY3(ntkl, cat, deg, k)==1) {
	ARRAY3(ch, cat, deg, k)=edges;
      }
    }
  }
  
  /* measurement done, update change */
  for (k=0; k<no_cats; k++) {
    for (i=0; i<maxind+1; i++) {
      for (j=0; j<agebins; j++) {
	igraph_real_t oldmean;
	if (ARRAY3(ntkl, k, i, j) !=0) {
	  ARRAY3(normfact, k, i, j) += (edges-ARRAY3(ch, k, i, j)+1);
	}
	oldmean=ARRAY3(*adkl, k, i, j);
	ARRAY3(*adkl, k, i, j) *=
	  ARRAY3(notnull, k, i, j)/ARRAY3(normfact, k, i, j);
	if (lsd) {
	  ARRAY3(*sd, k, i, j) += 
	    oldmean * oldmean * ARRAY3(notnull, k, i, j)*
	    (1-ARRAY3(notnull, k, i, j)/ARRAY3(normfact, k, i, j));
	  if (ARRAY3(normfact, k, i, j) > 0) {
	    ARRAY3(*sd, k, i, j)=
	      sqrt(ARRAY3(*sd, k, i, j)/(ARRAY3(normfact, k, i, j)-1));
	  }
	}
      }
    }
  }
  
  igraph_array3_destroy(&normfact);
  
  igraph_Free(indegree);
  igraph_array3_destroy(&ntkl);
  igraph_array3_destroy(&ch);
  igraph_array3_destroy(&notnull);
  igraph_vector_destroy(&neis);

  return 0;
}

int igraph_measure_dynamics_citedcat_id_age_st(const igraph_t *graph,
					       igraph_vector_t *res,
					       const igraph_array3_t *adkl,
					       const igraph_vector_t *cats, 
					       igraph_integer_t pno_cats) {
  long int agebins=igraph_array3_n(adkl, 3);
  long int no_of_nodes=igraph_vcount(graph);
  long int binwidth;
  
  int *indegree;
  igraph_vector_t neis;
  
  long int node;
  long int i, k;

  igraph_vector_init(&neis, 0);
  indegree=igraph_Calloc(no_of_nodes, int);
  binwidth=no_of_nodes/agebins+1;
  
  igraph_vector_resize(res, no_of_nodes);
  igraph_vector_null(res);
  VECTOR(*res)[0]=ARRAY3(*adkl, (long int)VECTOR(*cats)[0], 0, 0);
  
  for (node=1; node<no_of_nodes; node++) {
    long int cidx;
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    /* new node, aging */
    cidx=VECTOR(*cats)[node];
    VECTOR(*res)[node] = VECTOR(*res)[node-1] + ARRAY3(*adkl, cidx, 0, 0);
    for (k=1; node-binwidth*k+1 >= 1; k++) {
      long int shnode=node-binwidth*k;
      long int cat=VECTOR(*cats)[shnode];
      long int deg=indegree[shnode];
      VECTOR(*res)[node] += 
	-ARRAY3(*adkl, cat, deg, k-1)+ARRAY3(*adkl, cat, deg, k);
    }

    /* inspect the outgoing edges */
    igraph_neighbors(graph, &neis, node, IGRAPH_OUT);
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int cidx=VECTOR(*cats)[to];      
      long int xidx=indegree[to];
      long int yidx=(node-to)/binwidth;
      
      indegree[to]++;
      VECTOR(*res)[node] += 
	-ARRAY3(*adkl, cidx, xidx, yidx) + ARRAY3(*adkl, cidx, xidx+1, yidx);
    }
  }

  igraph_vector_destroy(&neis);
  igraph_Free(indegree);
  
  return 0;
}

int igraph_measure_dynamics_citingcat_id_age(const igraph_t *graph,
					     igraph_array3_t *adkl,
					     igraph_array3_t *sd,
					     const igraph_vector_t *st,
					     const igraph_vector_t *cats,
					     igraph_integer_t pno_cats,
					     igraph_integer_t pagebins,
					     igraph_integer_t pmaxind) {
  
  long int agebins=pagebins;
  long int maxind=pmaxind;
  long int no_cats=pno_cats;
  long int no_of_nodes=igraph_vcount(graph);
  long int binwidth;
  
  int *indegree;
  igraph_matrix_t ntkl;
  igraph_array3_t ch, normfact, notnull;
  igraph_vector_t edges;
  igraph_vector_t neis;
  
  long int node;
  long int i,j,k;
  
  igraph_bool_t lsd=(sd != 0);
  
  binwidth=no_of_nodes/agebins+1;
  
  igraph_vector_init(&neis, 0);
  indegree=igraph_Calloc(no_of_nodes, int);
  igraph_vector_init(&edges, no_cats);

  igraph_array3_resize(adkl, no_cats, maxind+1, agebins);
  igraph_array3_null(adkl);
  if (lsd) {
    igraph_array3_resize(sd, no_cats, maxind+1, agebins);
    igraph_array3_null(sd);
  }
  igraph_matrix_init(&ntkl, maxind+1, agebins);
  igraph_array3_init(&ch,       no_cats, maxind+1, agebins);
  igraph_array3_init(&normfact, no_cats, maxind+1, agebins);
  igraph_array3_init(&notnull,  no_cats, maxind+1, agebins);
  
  for (node=0; node < no_of_nodes; node++) {
    long int cidx=VECTOR(*cats)[node];
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    /* update A(g,k,l) */
    igraph_neighbors(graph, &neis, node, IGRAPH_OUT);
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=indegree[to];
      long int yidx=(node-to)/binwidth;
      
      double xk=VECTOR(*st)[node] / MATRIX(ntkl, xidx, yidx);
      double oldm=ARRAY3(*adkl, cidx, xidx, yidx);
      ARRAY3(notnull, cidx, xidx, yidx) += 1;
      ARRAY3(*adkl, cidx, xidx, yidx) += 
	(xk-oldm)/ARRAY3(notnull, cidx, xidx, yidx);
      if (lsd) {
	ARRAY3(*sd, cidx, xidx, yidx) += 
	  (xk-oldm)*(xk-ARRAY3(*adkl, cidx, xidx, yidx));
      }
    }
    
    /* add the new edges */
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=indegree[to];
      long int yidx=(node-to)/binwidth;
      
      indegree[to]++;
      MATRIX(ntkl, xidx, yidx)--;
      if (MATRIX(ntkl, xidx, yidx)==0) {
	for (j=0; j<no_cats; j++) {
	  ARRAY3(normfact, j, xidx, yidx) += 
	    (VECTOR(edges)[j]-ARRAY3(ch, j, xidx, yidx)+1);
	  ARRAY3(ch, j, xidx, yidx)=VECTOR(edges)[j];
	}
      }
      MATRIX(ntkl, xidx+1, yidx)++;
      if (MATRIX(ntkl, xidx+1, yidx)==1) {
	for (j=0; j<no_cats; j++) {
	  ARRAY3(ch, j, xidx+1, yidx)=VECTOR(edges)[j];
	}
      }
      VECTOR(edges)[cidx]++;
    }
    
    /* aging */
    MATRIX(ntkl, 0, 0)++;
    if (MATRIX(ntkl, 0, 0)==1) {
      for (j=0; j<no_cats; j++) {
	ARRAY3(ch, j, 0, 0)=VECTOR(edges)[j];
      }
    }
    for (k=1; node-binwidth*k+1 >= 1; k++) {
      long int shnode=node-binwidth*k;
      long int deg=indegree[shnode];
      MATRIX(ntkl, deg, k-1)--;
      if (MATRIX(ntkl, deg, k-1)==0) {
	for (j=0; j<no_cats; j++) {
	  ARRAY3(normfact, j, deg, k-1) += 
	    (VECTOR(edges)[j]-ARRAY3(ch, j, deg, k-1)+1);
	  ARRAY3(ch, j, deg, k-1)=VECTOR(edges)[j];
	}
      }
      MATRIX(ntkl, deg, k)++;
      if (MATRIX(ntkl, deg, k)==1) {
	for (j=0; j<no_cats; j++) {
	  ARRAY3(ch, j, deg, k)=VECTOR(edges)[j];
	}
      }
    }
  }
  
  /* measurement done, update change */
  for (k=0; k<no_cats; k++) {
    for (i=0; i<maxind+1; i++) {
      for (j=0; j<agebins; j++) {
	igraph_real_t oldmean;
	if (MATRIX(ntkl, i, j) != 0) {
	  ARRAY3(normfact, k, i, j) +=
	    (VECTOR(edges)[k]-ARRAY3(ch, k, i, j)+1);
	}
	oldmean=ARRAY3(*adkl, k, i, j);
	ARRAY3(*adkl, k, i, j) *= 
	  ARRAY3(notnull, k, i, j)/ARRAY3(normfact, k, i, j);
	if (lsd) {
	  ARRAY3(*sd, k, i, j) += 
	    oldmean * oldmean * ARRAY3(notnull, k, i, j)*
	    (1-ARRAY3(notnull, k, i, j)/ARRAY3(normfact, k, i, j));
	  if (ARRAY3(normfact, k, i, j) > 0) {
	    ARRAY3(*sd, k, i, j)=
	      sqrt(ARRAY3(*sd, k, i, j)/(ARRAY3(normfact, k, i, j)-1));
	  }
	}
      }
    }
  }
  
  igraph_array3_destroy(&normfact);
  igraph_Free(indegree);
  igraph_matrix_destroy(&ntkl);
  igraph_array3_destroy(&ch);
  igraph_array3_destroy(&notnull);
  igraph_vector_destroy(&neis);
  igraph_vector_destroy(&edges);
      
  return 0;
}

int igraph_measure_dynamics_citingcat_id_age_st(const igraph_t *graph,
						igraph_vector_t *res,
						const igraph_array3_t *adkl,
						const igraph_vector_t *cats,
						igraph_integer_t pno_cats) {
  
  long int agebins=igraph_array3_n(adkl, 3);
  long int no_of_nodes=igraph_vcount(graph);
  long int no_cats=pno_cats;
  long int binwidth;
  
  int *indegree;
  igraph_vector_t neis;
  
  long int node;
  long int i, j, k;
  
  igraph_matrix_t allst;
  
  igraph_matrix_init(&allst, no_cats, no_of_nodes+1);
  igraph_vector_init(&neis, 0);
  indegree=igraph_Calloc(no_of_nodes, int);
  binwidth=no_of_nodes/agebins+1;
  
  igraph_vector_resize(res, no_of_nodes);
  igraph_vector_null(res);
  for (j=0; j<no_cats; j++) {
    MATRIX(allst, j, 0)=ARRAY3(*adkl, j, 0, 0);
  }
  VECTOR(*res)[0]=MATRIX(allst, (long int)VECTOR(*cats)[0], 0);
  
  for (node=1; node<no_of_nodes; node++) {
    long int cidx=VECTOR(*cats)[node];
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    /* new node, aging */
    for (j=0; j<no_cats; j++) {
      MATRIX(allst, j, node) = MATRIX(allst, j, node-1) + 
	ARRAY3(*adkl, j, 0, 0);
    }
    for (k=1; node-binwidth*k+1 >= 1; k++) {
      long int shnode=node-binwidth*k;
      long int deg=indegree[shnode];
      for (j=0; j<no_cats; j++) {
	MATRIX(allst, j, node) += 
	  -ARRAY3(*adkl, j, deg, k-1) + ARRAY3(*adkl, j, deg, k);
      }
    }

    /* inspect the outgoing edges */
    igraph_neighbors(graph, &neis, node, IGRAPH_OUT);
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=indegree[to];
      long int yidx=(node-to)/binwidth;
      
      indegree[to]++;
      for (j=0; j<no_cats; j++) {
	MATRIX(allst, j, node) +=
	  -ARRAY3(*adkl, j, xidx, yidx) + ARRAY3(*adkl, j, xidx+1, yidx);
      }
    }

    /* update the result */
    VECTOR(*res)[node]=MATRIX(allst, cidx, node);
  }

  igraph_vector_destroy(&neis);
  igraph_matrix_destroy(&allst);
  igraph_Free(indegree);
  
  return 0;
}

#define NTKK(xidx, yidx) \
   ((xidx)==(yidx)) ? (VECTOR(ntk)[(xidx)]*(VECTOR(ntk)[(xidx)]-1)/2-MATRIX(ntkk,(xidx),(yidx))) : (VECTOR(ntk)[(xidx)]*VECTOR(ntk)[(yidx)]-MATRIX(ntkk,(xidx),(yidx)))

/* int print_ntkk(igraph_matrix_t *ntkk, igraph_vector_t *ntk) { */
/*   long int i, j, r=igraph_matrix_nrow(ntkk), c=igraph_matrix_ncol(ntkk); */
/*   for (i=0; i<r; i++) { */
/*     for (j=0; j<c; j++) { */
/*       long int val=(i==j) ?  */
/* 	(VECTOR(*ntk)[i]*(VECTOR(*ntk)[i]-1)/2-MATRIX(*ntkk,i,j)) :  */
/* 	(VECTOR(*ntk)[i]*VECTOR(*ntk)[j]-MATRIX(*ntkk,i,j)); */
/*       fprintf(stderr, "%li ", val); */
/*     } */
/*     fprintf(stderr, "\n"); */
/*   } */
/*   return 0; */
/* } */

int igraph_measure_dynamics_d_d(const igraph_t *graph,
				const igraph_vector_t *ntime,
				const igraph_vector_t *etime,
				igraph_integer_t events,
				igraph_matrix_t *akk,
				igraph_matrix_t *sd,
				const igraph_vector_t *st,
				igraph_integer_t pmaxdeg) {
  
  long int maxdeg=pmaxdeg;
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);

  igraph_vector_t adjedges;
  
  igraph_vector_t degree;	/* the actual degrees of the nodes */
  igraph_vector_t ntk;		/* the actual number of k-nodes */
  igraph_matrix_t ntkk;         /* the number of (k1-k2) edges */
  igraph_matrix_t normfact;	/* number of experiments */
  igraph_matrix_t ch;		/* the time of the last switch */
  igraph_matrix_t notnull;	/* the number of non-zero experiments */
  igraph_vector_t added;	/* whether the edge is added already */

  igraph_vector_t ntimeidx;	/* lookup vector for nodes */
  igraph_vector_t etimeidx;     /* lookup vector for edges */

  long int timestep=0;
  long int nptr=0;
  long int eptr=0;
  long int eptr_save, nptr_save, i, j;

  /* Init everything */
  IGRAPH_VECTOR_INIT_FINALLY(&adjedges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&ntk, maxdeg+1);
  IGRAPH_MATRIX_INIT_FINALLY(&ntkk, maxdeg+1, maxdeg+1);
  IGRAPH_MATRIX_INIT_FINALLY(&normfact, maxdeg+1, maxdeg+1);
  IGRAPH_MATRIX_INIT_FINALLY(&ch, maxdeg+1, maxdeg+1);
  IGRAPH_MATRIX_INIT_FINALLY(&notnull, maxdeg+1, maxdeg+1);
  IGRAPH_VECTOR_INIT_FINALLY(&added, no_of_edges);
  IGRAPH_VECTOR_INIT_FINALLY(&ntimeidx, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&etimeidx, 0);

  /* Resize and init the output */
  IGRAPH_CHECK(igraph_matrix_resize(akk, maxdeg+1, maxdeg+1));
  igraph_matrix_null(akk);
  if (sd) {
    IGRAPH_CHECK(igraph_matrix_resize(sd, maxdeg+1, maxdeg+1));
    igraph_matrix_null(sd);
  }
  
  /* Create lookup vectors for nodes and edges */
  IGRAPH_CHECK(igraph_vector_order1(ntime, &ntimeidx, events));
  IGRAPH_CHECK(igraph_vector_order1(etime, &etimeidx, events));

  for (timestep=0; timestep<events; timestep++) {

    IGRAPH_ALLOW_INTERRUPTION();

/*     fprintf(stderr, "-----------time step %li\n", timestep); */

    /* Add the nodes */

    nptr_save=nptr;
    while (nptr < no_of_nodes &&
	   VECTOR(*ntime)[ (long int) VECTOR(ntimeidx)[nptr] ] == timestep) {
      nptr++;
    }
    VECTOR(ntk)[0] += (nptr-nptr_save);
    if (VECTOR(ntk)[0] == nptr-nptr_save && nptr-nptr_save > 0) {
      /* just introduced zero degree nodes to the net, update ch */
      for (i=0; i<maxdeg+1; i++) {
	MATRIX(ch, 0, i) = MATRIX(ch, i, 0) = eptr;
      }
    }
    
/*     print_ntkk(&ntkk, &ntk); */
    for (i=0; i<maxdeg+1; i++) {
      for (j=0; j<=i; j++) {
	if (NTKK(i,j) > 0) {
	  MATRIX(normfact, i, j) += 1;
	  MATRIX(normfact, j, i) = MATRIX(normfact, i, j);
	}
      }
    }
    
    /* Estimage Akk */

    eptr_save=eptr;
    while (eptr < no_of_edges &&
	   VECTOR(*etime)[ (long int) VECTOR(etimeidx)[eptr] ] == timestep) {
      long int edge=VECTOR(etimeidx)[eptr];
      long int xidx, yidx;
      igraph_integer_t from, to;
      double xk, oldakk, ntkkval;

      igraph_edge(graph, edge, &from, &to);
      xidx=VECTOR(degree)[ (long int)from ];
      yidx=VECTOR(degree)[ (long int)to ];
      MATRIX(notnull, xidx, yidx) += 1;
      MATRIX(notnull, yidx, xidx) += 1;	/* bug */

      ntkkval=NTKK(xidx, yidx);
      xk=VECTOR(*st)[timestep]/ntkkval;
      oldakk=MATRIX(*akk, xidx, yidx);
      MATRIX(*akk, xidx, yidx) = oldakk +
	(xk-oldakk)/(MATRIX(notnull, xidx, yidx));
      MATRIX(*akk, yidx, xidx) = MATRIX(*akk, xidx, yidx);
      if (sd) {
	MATRIX(*sd, xidx, yidx) += (xk-oldakk)*(xk-MATRIX(*akk, xidx, yidx));
	MATRIX(*sd, yidx, xidx) = MATRIX(*sd, xidx, yidx);
      }
      eptr++;
    }

    /* Add the edges */
    
    eptr=eptr_save;
    while (eptr < no_of_edges &&
	   VECTOR(*etime)[ (long int) VECTOR(etimeidx)[eptr] ] == timestep) {
      long int edge=VECTOR(etimeidx)[eptr];
      long int xidx, yidx;
      igraph_integer_t from, to;
      
      igraph_edge(graph, edge, &from, &to);
      xidx=VECTOR(degree)[ (long int)from ];
      yidx=VECTOR(degree)[ (long int)to ];

      igraph_adjacent(graph, &adjedges, from, IGRAPH_ALL);
      for (i=0; i<igraph_vector_size(&adjedges); i++) {
	igraph_integer_t e_from, e_to;
	long int deg;
	long int edge=VECTOR(adjedges)[i];
	igraph_edge(graph, edge, &e_from, &e_to);
	if (VECTOR(added)[edge]) {
	  if (e_to==from) { e_to=e_from; }
	  deg=VECTOR(degree)[(long int)e_to];
	  MATRIX(ntkk, xidx, deg) =MATRIX(ntkk, xidx, deg)-1;
	  MATRIX(ntkk, deg, xidx) =MATRIX(ntkk, xidx, deg);
	  MATRIX(ntkk, xidx+1, deg) =MATRIX(ntkk, xidx+1, deg)+1;
	  MATRIX(ntkk, deg, xidx+1) =MATRIX(ntkk, xidx+1, deg);
	}
      }
      igraph_adjacent(graph, &adjedges, to, IGRAPH_ALL);
      for (i=0; i<igraph_vector_size(&adjedges); i++) {
	igraph_integer_t e_from, e_to;
	long int deg;
	long int edge=VECTOR(adjedges)[i];
	igraph_edge(graph, edge, &e_from, &e_to);
	if (VECTOR(added)[edge]) {
	  if (e_to==to) { e_to=e_from; }
	  deg=VECTOR(degree)[(long int)e_to];
	  MATRIX(ntkk, yidx, deg) = MATRIX(ntkk, yidx, deg) - 1;
	  MATRIX(ntkk, deg, yidx) = MATRIX(ntkk, yidx, deg);
	  MATRIX(ntkk, yidx+1, deg) =MATRIX(ntkk, yidx+1, deg)+1;
	  MATRIX(ntkk, deg, yidx+1) = MATRIX(ntkk, yidx+1, deg);
	}
      }
      MATRIX(ntkk, xidx+1, yidx+1) += 1;
      MATRIX(ntkk, yidx+1, xidx+1) =MATRIX(ntkk, xidx+1, yidx+1);
      VECTOR(added)[edge]=1;
      
      VECTOR(ntk)[xidx] --;
      VECTOR(ntk)[yidx] --;
      VECTOR(ntk)[xidx+1] ++;
      VECTOR(ntk)[yidx+1] ++;
      
      VECTOR(degree)[ (long int)from ] ++;
      VECTOR(degree)[ (long int)to ] ++;
      
      eptr++;
    }

/*     fprintf(stderr, "--------\n"); */
/*     print_ntkk(&ntkk, &ntk); */

  }

/*   for (i=0; i<maxdeg+1; i++) { */
/*     if (VECTOR(ntk)[i] != 0) { */
/*       for (j=0; j<=i; j++) { */
/* 	if (NTKK(i, j) > 0) { */
/* 	  MATRIX(normfact, i, j) += eptr-MATRIX(ch, i, j); */
/* 	  MATRIX(normfact, j, i) = MATRIX(normfact, i, j); */
/* 	} */
/*       } */
/*     } */
/*   } */
  
/*   fprintf(stderr, "---------\n"); */
/*   print_matrix(&normfact); */

  /* Update akk, sd */
  for (i=0; i<maxdeg+1; i++) {
    igraph_real_t oldakk;
    for (j=0; j<=i; j++) {
      oldakk=MATRIX(*akk, i, j);
      MATRIX(*akk, i, j) *= MATRIX(notnull, i, j) / MATRIX(normfact, i, j);
      MATRIX(*akk, j, i) = MATRIX(*akk, i, j);
      if (sd) {
	MATRIX(*sd, i, j) += oldakk * oldakk * MATRIX(notnull, i, j) *
	  (1-MATRIX(notnull, i, j)/MATRIX(normfact, i, j));
	if (MATRIX(normfact, i, j) > 0) {
	  MATRIX(*sd, i, j)=sqrt(MATRIX(*sd, i, j)/(MATRIX(normfact, i, j)-1));
	  MATRIX(*sd, j, i) = MATRIX(*sd, i, j);
	}
      }
    }
  }

  igraph_vector_destroy(&etimeidx);
  igraph_vector_destroy(&ntimeidx);
  igraph_vector_destroy(&added);
  igraph_matrix_destroy(&notnull);
  igraph_matrix_destroy(&ch);
  igraph_matrix_destroy(&normfact);
  igraph_matrix_destroy(&ntkk);
  igraph_vector_destroy(&ntk);
  igraph_vector_destroy(&degree);
  igraph_vector_destroy(&adjedges);

  IGRAPH_FINALLY_CLEAN(10);
  return 0;
}  
				  
int igraph_measure_dynamics_d_d_st(const igraph_t *graph,        /* input */
				   const igraph_vector_t *ntime, /* input */
				   const igraph_vector_t *etime, /* input */
				   const igraph_matrix_t *akk,   /* input */
				   igraph_integer_t events,
				   igraph_integer_t maxtotaldeg,
				   igraph_vector_t *st) {        /* output */
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  long int timestep=0;
  long int nptr=0;
  long int eptr=0;
  long int maxdeg=0;

  igraph_vector_t degree;
  igraph_vector_t ntk;		/* number of nodes of different types */
  igraph_vector_t added;
  igraph_vector_t adjedges;
  igraph_vector_t ntimeidx;	/* lookup vector for nodes */
  igraph_vector_t etimeidx;	/* lookup vector for edges */

  long int i;

  /* Init everything */
  IGRAPH_VECTOR_INIT_FINALLY(&ntk, maxtotaldeg+1);
  IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&added, no_of_edges);
  IGRAPH_VECTOR_INIT_FINALLY(&adjedges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&ntimeidx, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&etimeidx, 0);
  
  /* Resize result */
  IGRAPH_CHECK(igraph_vector_resize(st, events+1));
  VECTOR(*st)[0]=0;
  
  /* Create lookup vectors */
  IGRAPH_CHECK(igraph_vector_order1(ntime, &ntimeidx, events));
  IGRAPH_CHECK(igraph_vector_order1(etime, &etimeidx, events));

  for (timestep=0; timestep<events; timestep++) {

    IGRAPH_ALLOW_INTERRUPTION();
    
    /* add the new nodes, if any */
    while (nptr < no_of_nodes &&
	   VECTOR(*ntime)[ (long int) VECTOR(ntimeidx)[nptr] ] == timestep) {
      igraph_real_t akk_inc;
      akk_inc=0;
      for (i=0; i<=maxdeg; i++) {
	akk_inc += VECTOR(ntk)[i]*MATRIX(*akk, i, 0);
	akk_inc += VECTOR(ntk)[i]*MATRIX(*akk, i, 0); /* why this? */
      }
      VECTOR(*st)[timestep] += akk_inc;
      VECTOR(ntk)[0]++;

      nptr++;
    }

    /* add the new edges if any */
    while (eptr < no_of_edges &&
	   VECTOR(*etime)[ (long int) VECTOR(etimeidx)[eptr] ] == timestep) {
      igraph_real_t akk_inc;
      long int edge=VECTOR(etimeidx)[eptr];
      igraph_integer_t from, to;
      long int xidx, yidx;
      igraph_edge(graph, edge, &from, &to);
      xidx=VECTOR(degree)[ (long int)from];
      yidx=VECTOR(degree)[ (long int)to];

      VECTOR(degree)[(long int)from]++;
      VECTOR(degree)[(long int)to]++;

      akk_inc=0;

      for (i=0; i<=maxdeg+1; i++) {
	akk_inc -= VECTOR(ntk)[i] * MATRIX(*akk, i, xidx);
	akk_inc -= VECTOR(ntk)[i] * MATRIX(*akk, i, yidx);
	akk_inc += VECTOR(ntk)[i] * MATRIX(*akk, i, xidx+1);
	akk_inc += VECTOR(ntk)[i] * MATRIX(*akk, i, yidx+1);
      }
      akk_inc += MATRIX(*akk, xidx, xidx);
      akk_inc += MATRIX(*akk, yidx, yidx); /* why twice???? */
      akk_inc -= MATRIX(*akk, xidx+1, xidx+1);
      akk_inc -= MATRIX(*akk, yidx+1, yidx+1);

      VECTOR(ntk)[xidx]--;
      VECTOR(ntk)[yidx]--;
      VECTOR(ntk)[xidx+1]++;
      VECTOR(ntk)[yidx+1]++;
      for (i=0; i<=maxdeg; i++) {
	akk_inc += VECTOR(ntk)[i]*
	  (MATRIX(*akk, i, xidx+1)-MATRIX(*akk, i, xidx)+
	   MATRIX(*akk, i, yidx+1)-MATRIX(*akk, i, yidx));
      }

      if (xidx+1 > maxdeg) {
	maxdeg=xidx+1;
      }
      if (yidx+1 > maxdeg) {
	maxdeg=yidx+1;
      }

      /* Now update the edges */
      igraph_adjacent(graph, &adjedges, from, IGRAPH_ALL);
      for (i=0; i<igraph_vector_size(&adjedges); i++) {
	igraph_integer_t e_from, e_to;
	long int deg;
	long int edge=VECTOR(adjedges)[i];
	igraph_edge(graph, edge, &e_from, &e_to);
	if (VECTOR(added)[edge]) {
	  if (e_to==from) { e_to=e_from; }
	  deg=VECTOR(degree)[(long int)e_to];
	  akk_inc += MATRIX(*akk, xidx, deg);
	  akk_inc -= MATRIX(*akk, xidx+1, deg);
	}
      }
      igraph_adjacent(graph, &adjedges, to, IGRAPH_ALL);
      for (i=0; i<igraph_vector_size(&adjedges); i++) {
	igraph_integer_t e_from, e_to;
	long int deg;
	long int edge=VECTOR(adjedges)[i];
	igraph_edge(graph, edge, &e_from, &e_to);
	if (VECTOR(added)[edge]) {
	  if (e_to==to) { e_to=e_from; }
	  deg=VECTOR(degree)[(long int)e_to];
	  akk_inc += MATRIX(*akk, yidx, deg);
	  akk_inc -= MATRIX(*akk, yidx+1, deg);
	}
      }
      VECTOR(added)[edge]=1;

      VECTOR(*st)[timestep] += akk_inc;
      eptr++;
    }
    
    VECTOR(*st)[timestep+1]=VECTOR(*st)[timestep];       
  }

  igraph_vector_pop_back(st);
  
  igraph_vector_destroy(&etimeidx);
  igraph_vector_destroy(&ntimeidx);
  igraph_vector_destroy(&adjedges);
  igraph_vector_destroy(&added);
  igraph_vector_destroy(&degree);
  igraph_vector_destroy(&ntk);
  IGRAPH_FINALLY_CLEAN(6);
  
  return 0;
}

/* int print_vector(const igraph_vector_t *v) { */
/*   long int i, n=igraph_vector_size(v); */
/*   for (i=0; i<n; i++) { */
/*     fprintf(stderr, "%li ", (long int) VECTOR(*v)[i]); */
/*   } */
/*   fprintf(stderr, "\n"); */
/*   return 0; */
/* } */

/* int print_matrix(const igraph_matrix_t *m) { */
/*   long int i, j, r=igraph_matrix_nrow(m), c=igraph_matrix_ncol(m); */
/*   for (i=0; i<r; i++) { */
/*     for (j=0; j<c; j++) { */
/*       fprintf(stderr, "%li ", (long int) MATRIX(*m, i, j)); */
/*     } */
/*     fprintf(stderr, "\n"); */
/*   } */
/*   return 0; */
/* } */

int igraph_measure_dynamics_lastcit(const igraph_t *graph, igraph_vector_t *al,
				    igraph_vector_t *sd, 
				    igraph_vector_t *no,
				    const igraph_vector_t *st,
				    igraph_integer_t pagebins) {

  long int no_of_nodes=igraph_vcount(graph);

  /* there will be 'agebins' categories, plus one for vertices without a 
     citation.
  */
  long int agebins=pagebins;
  long int binwidth=no_of_nodes/agebins+1;

  igraph_vector_t ntl, ch, notnull, normfact;

  /* lastcit[n] is the time step+1 in which n was cited, if not 0;
     if 0, that means that n was never cited
   */
  long int *lastcit;

  long int node, i, k;

  igraph_vector_t neis;
  
  long int edges=0;
  
  lastcit=igraph_Calloc(no_of_nodes, long int);
  if (!lastcit) {
    IGRAPH_ERROR("Cannot measure dynamics (lastcit)", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, lastcit);

  IGRAPH_VECTOR_INIT_FINALLY(&ntl, agebins+1);
  IGRAPH_VECTOR_INIT_FINALLY(&ch, agebins+1);
  IGRAPH_VECTOR_INIT_FINALLY(&notnull, agebins+1);
  IGRAPH_VECTOR_INIT_FINALLY(&normfact, agebins+1);

  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  
  IGRAPH_CHECK(igraph_vector_resize(al, agebins+1));
  igraph_vector_null(al);
  if (sd) {
    IGRAPH_CHECK(igraph_vector_resize(sd, agebins+1));
    igraph_vector_null(sd);
  }
  
  for (node=0; node<no_of_nodes; node++) {

    IGRAPH_ALLOW_INTERRUPTION();

    /* update A(.) */
    igraph_neighbors(graph, &neis, node, IGRAPH_OUT);
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=lastcit[to]!=0 ? (node-lastcit[to]+1)/binwidth : agebins;
      
      double xk=VECTOR(*st)[node]/VECTOR(ntl)[xidx];
      double oldm=VECTOR(*al)[xidx];

      VECTOR(notnull)[xidx] += 1;
      VECTOR(*al)[xidx] += (xk-oldm)/VECTOR(notnull)[xidx];
      if (sd) {
	VECTOR(*sd)[xidx] += (xk-oldm)*(xk-VECTOR(*al)[xidx]);
      }
    }

    /* Update ntk, ch, normfact */
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=lastcit[to]!=0 ? (node-lastcit[to]+1)/binwidth : agebins;

      lastcit[to]=node+1;
      VECTOR(ntl)[xidx]--;
      if (VECTOR(ntl)[xidx]==0) {
	VECTOR(normfact)[xidx] += (edges-VECTOR(ch)[xidx]+1);
	VECTOR(ch)[xidx]=edges;
      }
      VECTOR(ntl)[0]++;
      if (VECTOR(ntl)[0]==1) {
	VECTOR(ch)[0]=edges;
      }
      edges++;
    }    

    /* the new node */
    VECTOR(ntl)[agebins]++;
    if (VECTOR(ntl)[agebins]==1) {
      VECTOR(ch)[agebins]=edges;
    }
    
    /* should we move some citations to an older bin? */
    for (k=1; node-binwidth*k+1 >= 1; k++) {
      long int shnode=node-binwidth*k;
      igraph_neighbors(graph, &neis, shnode, IGRAPH_OUT);
      for (i=0; i<igraph_vector_size(&neis); i++) {
	long int cnode=VECTOR(neis)[i];
	if (lastcit[cnode]==shnode+1) {
	  /* the last citation to cnode was made by shnode, move to other bin*/
	  VECTOR(ntl)[k-1]--;
	  if (VECTOR(ntl)[k-1]==0) {
	    VECTOR(normfact)[k-1] += (edges-VECTOR(ch)[k-1]+1);
	    VECTOR(ch)[k-1] = edges;
	  }
	  VECTOR(ntl)[k]++;
	  if (VECTOR(ntl)[k]==1) {
	    VECTOR(ch)[k]=edges;
	  }
	}
      }
    }
    
  }

  /* measurement done, update change */
  for (i=0; i<agebins+1; i++) {
    igraph_real_t oldmean;
    if (VECTOR(ntl)[i] != 0) {
      VECTOR(normfact)[i] += (edges-VECTOR(ch)[i]+1);
    }
    oldmean=VECTOR(*al)[i];
    VECTOR(*al)[i] *= VECTOR(notnull)[i]/VECTOR(normfact)[i];
    if (sd) {
      VECTOR(*sd)[i] += oldmean * oldmean * VECTOR(notnull)[i] *
	(1-VECTOR(notnull)[i]/VECTOR(normfact)[i]);
      if (VECTOR(normfact)[i] > 0) {
	VECTOR(*sd)[i] = sqrt(VECTOR(*sd)[i]/(VECTOR(normfact)[i]-1));
      }
    }
  }  

  if (no) {
    igraph_vector_destroy(no);
    *no=normfact;
  } else {
    igraph_vector_destroy(&normfact);
  }
 
  igraph_free(lastcit);
  igraph_vector_destroy(&ntl);
  igraph_vector_destroy(&ch);
  igraph_vector_destroy(&notnull);
  igraph_vector_destroy(&neis);
  IGRAPH_FINALLY_CLEAN(6);
  return 0;
}

int igraph_measure_dynamics_lastcit_st(const igraph_t *graph, 
				       igraph_vector_t *res,
				       const igraph_vector_t *al) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int agebins=igraph_vector_size(al)-1;
  long int binwidth=no_of_nodes/agebins+1;
  long int *lastcit;

  igraph_vector_t neis;

  long int node, i, k;

  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  lastcit=igraph_Calloc(no_of_nodes, long int);
  if (!lastcit) {
    IGRAPH_ERROR("Cannot measure dynamics (lastcit st)", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, lastcit);
  
  igraph_vector_resize(res, no_of_nodes);
  igraph_vector_null(res);
  VECTOR(*res)[0]=VECTOR(*al)[agebins];	/* node without citation */
  
  for (node=1; node<no_of_nodes; node++) {
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    /* new node, updates */
    VECTOR(*res)[node]=VECTOR(*res)[node-1]+VECTOR(*al)[agebins];
    for (k=1; node-binwidth*k+1 >= 1; k++) {
      long int shnode=node-binwidth*k;
      igraph_neighbors(graph, &neis, shnode, IGRAPH_OUT);
      for (i=0; i<igraph_vector_size(&neis); i++) {
	long int cnode=VECTOR(neis)[i];
	if (lastcit[cnode]==shnode+1) {
	  VECTOR(*res)[node] += -VECTOR(*al)[k-1]+VECTOR(*al)[k];
	}
      }
    }
    
    /* inspect the outgoing edges */
    igraph_neighbors(graph, &neis, node, IGRAPH_OUT);
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=lastcit[to]!=0 ? (node-lastcit[to]+1)/binwidth : agebins;
      
      lastcit[to]=node+1;
      
      VECTOR(*res)[node] += -VECTOR(*al)[xidx]+VECTOR(*al)[0];
    }
  }
  
  igraph_free(lastcit);
  igraph_vector_destroy(&neis);
  IGRAPH_FINALLY_CLEAN(2);

  return 0;
}

int igraph_measure_dynamics_age(const igraph_t *graph, 
				igraph_vector_t *al,
				igraph_vector_t *sd,
				igraph_vector_t *no,
				const igraph_vector_t *st,
				igraph_integer_t pagebins) {

  long int no_of_nodes=igraph_vcount(graph);
  
  long int agebins=pagebins;
  long int binwidth=no_of_nodes/agebins+1;

  igraph_vector_t ntl, ch, notnull, normfact;
  
  long int node, i, k;
  
  igraph_vector_t neis;
  
  long int edges=0;
  
  IGRAPH_VECTOR_INIT_FINALLY(&ntl, agebins);
  IGRAPH_VECTOR_INIT_FINALLY(&ch, agebins);
  IGRAPH_VECTOR_INIT_FINALLY(&notnull, agebins);
  IGRAPH_VECTOR_INIT_FINALLY(&normfact, agebins);
  
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  
  IGRAPH_CHECK(igraph_vector_resize(al, agebins));
  igraph_vector_null(al);
  if (sd) {
    IGRAPH_CHECK(igraph_vector_resize(sd, agebins));
    igraph_vector_null(sd);
  }
  
  for (node=0; node<no_of_nodes; node++) {
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    /* update A() */
    igraph_neighbors(graph, &neis, node, IGRAPH_OUT);
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=(node-to)/binwidth;
      
      double xk=VECTOR(*st)[node]/VECTOR(ntl)[xidx];
      double oldm=VECTOR(*al)[xidx];
      
      VECTOR(notnull)[xidx] += 1;
      VECTOR(*al)[xidx] += (xk-oldm)/VECTOR(notnull)[xidx];
      if (sd) {
	VECTOR(*sd)[xidx] += (xk-oldm)*(xk-VECTOR(*al)[xidx]);
      }
    }
    
    /* update ntl, ch, normfact */
    /* these is omitted, no change in ntl, except for aging */
    edges += igraph_vector_size(&neis);
    
    /* new node */
    VECTOR(ntl)[0]++;
    if (VECTOR(ntl)[0]==1) {
      VECTOR(ch)[0]=edges;
    }
    
    /* aging, this could be written much simpler */
    for (k=1; node-binwidth*k+1 >=1; k++) {
      VECTOR(ntl)[k-1]--;
      if (VECTOR(ntl)[k-1]==0) {
	VECTOR(normfact)[k-1] += (edges-VECTOR(ch)[k-1]+1);
	VECTOR(ch)[k-1] = edges;
      }
      VECTOR(ntl)[k]++;
      if (VECTOR(ntl)[k]==1) {
	VECTOR(ch)[k]=edges;
      }
    }
    
  } /* node<no_of_nodes */

  /* measurement done, update change */
  for (i=0; i<agebins; i++) {
    igraph_real_t oldmean;
    if (VECTOR(ntl)[i] != 0) {
      VECTOR(normfact)[i] += (edges-VECTOR(ch)[i]+1);
    }
    oldmean=VECTOR(*al)[i];
    VECTOR(*al)[i] *= VECTOR(notnull)[i]/VECTOR(normfact)[i];
    if (sd) {
      VECTOR(*sd)[i] += oldmean * oldmean * VECTOR(notnull)[i] *
	(1-VECTOR(notnull)[i]/VECTOR(normfact)[i]);
      if (VECTOR(normfact)[i] > 0) {
	VECTOR(*sd)[i] = sqrt(VECTOR(*sd)[i]/(VECTOR(normfact)[i]-1));
      }
    }
  }
  
  if (no) {
    igraph_vector_destroy(no);
    *no=normfact;
  } else {
    igraph_vector_destroy(&normfact);
  }

  igraph_vector_destroy(&ntl);
  igraph_vector_destroy(&ch);
  igraph_vector_destroy(&notnull);
  igraph_vector_destroy(&neis);
  IGRAPH_FINALLY_CLEAN(5);
  return 0;
}

int igraph_measure_dynamics_age_st(const igraph_t *graph, 
				   igraph_vector_t *res,
				   const igraph_vector_t *al) {
  
  long int agebins=igraph_vector_size(al);
  long int no_of_nodes=igraph_vcount(graph);
  long int binwidth=no_of_nodes/agebins+1;
  
  long int node, k;
  
  IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
  igraph_vector_null(res);
  VECTOR(*res)[0]=VECTOR(*al)[0];
  
  for (node=1; node<no_of_nodes; node++) {
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    /* new node, aging */
    VECTOR(*res)[node] = VECTOR(*res)[node-1] + VECTOR(*al)[0];
    for (k=1; node-binwidth*k+1 >= 1; k++) {
      VECTOR(*res)[node] += -VECTOR(*al)[k-1] + VECTOR(*al)[k];
    }
    
    /* inspecting outgoing edges is not needed at all */

  }
  
  return 0;
}

int igraph_measure_dynamics_citedcat(const igraph_t *graph, 
				     const igraph_vector_t *cats,
				     igraph_integer_t pnocats,
				     igraph_vector_t *ak, 
				     igraph_vector_t *sd,
				     igraph_vector_t *no,
				     const igraph_vector_t *st) {
  
  long int nocats=pnocats;
  long int no_of_nodes=igraph_vcount(graph);
  
  igraph_vector_t normfact, ntk, ch, notnull;
  igraph_vector_t neis;
  
  long int node, i;
  long int edges=0;
  
  IGRAPH_VECTOR_INIT_FINALLY(&normfact, nocats);
  IGRAPH_VECTOR_INIT_FINALLY(&ntk, nocats);
  IGRAPH_VECTOR_INIT_FINALLY(&ch, nocats);
  IGRAPH_VECTOR_INIT_FINALLY(&notnull, nocats);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);

  IGRAPH_CHECK(igraph_vector_resize(ak, nocats));
  igraph_vector_null(ak);
  if (sd) {
    IGRAPH_CHECK(igraph_vector_resize(sd, nocats));
    igraph_vector_null(sd);
  }
  
  for (node=0; node<no_of_nodes; node++) {
    
    long int ccat=VECTOR(*cats)[node];
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    /* estimake A */
    igraph_neighbors(graph, &neis, node, IGRAPH_OUT);
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=VECTOR(*cats)[to];
      
      double xk=VECTOR(*st)[node]/VECTOR(ntk)[xidx];
      double oldm=VECTOR(*ak)[xidx];
      VECTOR(notnull)[xidx] += 1;
      VECTOR(*ak)[xidx] += (xk-oldm)/VECTOR(notnull)[xidx];
      if (sd) {
	VECTOR(*sd)[xidx] += (xk-oldm)*(xk-VECTOR(*ak)[xidx]);
      }
    }

    /* update ntk, ch, normfact, only the new node is important */
    edges += igraph_vector_size(&neis);
        
    VECTOR(ntk)[ccat] += 1;
    if (VECTOR(ntk)[ccat] == 1) {
      VECTOR(ch)[ccat]=edges;
    }
    
  }

  /* measurement done, update ch, normfact and ak and sd */
  for (i=0; i<nocats; i++) {
    igraph_real_t oldmean;
    if (VECTOR(ntk)[i] != 0) {
      VECTOR(normfact)[i] += (edges-VECTOR(ch)[i]+1);
    }
    oldmean=VECTOR(*ak)[i];
    VECTOR(*ak)[i] *= VECTOR(notnull)[i] / VECTOR(normfact)[i];
    if (sd) {
      VECTOR(*sd)[i] += oldmean * oldmean * VECTOR(notnull)[i] *
	(1-VECTOR(notnull)[i]/VECTOR(normfact)[i]);
      if (VECTOR(normfact)[i]>0) {
	VECTOR(*sd)[i] = sqrt(VECTOR(*sd)[i] / (VECTOR(normfact)[i]-1));
      }
    }
  }
  
  if (no) {
    igraph_vector_destroy(no);
    *no=normfact;
  } else {
    igraph_vector_destroy(&normfact);
  }  

  igraph_vector_destroy(&neis);
  igraph_vector_destroy(&ntk);
  igraph_vector_destroy(&ch);
  igraph_vector_destroy(&notnull);
  IGRAPH_FINALLY_CLEAN(5);
  
  return 0;
}

int igraph_measure_dynamics_citedcat_st(const igraph_t *graph,
					igraph_vector_t *res,
					const igraph_vector_t *ak,
					const igraph_vector_t *cats,
					igraph_integer_t pnocats) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int node;

  IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
  VECTOR(*res)[0]=VECTOR(*ak)[ (long int) VECTOR(*cats)[0] ];
  
  for (node=1; node<no_of_nodes; node++) {
    long int cidx=VECTOR(*cats)[node];
    IGRAPH_ALLOW_INTERRUPTION();
    VECTOR(*res)[node]=VECTOR(*res)[node-1] + VECTOR(*ak)[cidx];
  }    
  
  return 0;
}

int igraph_measure_dynamics_citingcat_citedcat(const igraph_t *graph,
					       igraph_matrix_t *agd,
					       igraph_matrix_t *sd,
					       igraph_matrix_t *no,
					       const igraph_vector_t *st,
					       const igraph_vector_t *cats,
					       igraph_integer_t pnocats) {
  long int nocats=pnocats;
  long int no_of_nodes=igraph_vcount(graph);
  
  igraph_vector_t ntd;
  igraph_matrix_t ch, normfact, notnull;
  igraph_vector_t neis;
  
  long int node;
  long int i, j;
  igraph_vector_t edges;
  
  IGRAPH_VECTOR_INIT_FINALLY(&ntd, nocats);
  IGRAPH_MATRIX_INIT_FINALLY(&ch, nocats, nocats);
  IGRAPH_MATRIX_INIT_FINALLY(&normfact, nocats, nocats);
  IGRAPH_MATRIX_INIT_FINALLY(&notnull, nocats, nocats);  
  
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&edges, nocats);
  
  IGRAPH_CHECK(igraph_matrix_resize(agd, nocats, nocats));
  igraph_matrix_null(agd);
  if (sd) {
    IGRAPH_CHECK(igraph_matrix_resize(sd, nocats, nocats));
    igraph_matrix_null(sd);
  }
  
  for (node=0; node<no_of_nodes; node++) {
    long int citingcat=VECTOR(*cats)[node];
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    /* update A() */
    igraph_neighbors(graph, &neis, node, IGRAPH_OUT);
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=VECTOR(*cats)[to];
      
      double xk=VECTOR(*st)[node] / VECTOR(ntd)[xidx];
      double oldm=MATRIX(*agd, citingcat, xidx);
      MATRIX(notnull, citingcat, xidx) += 1;
      MATRIX(*agd, citingcat, xidx) += 
	(xk-oldm)/MATRIX(notnull, citingcat, xidx);
      if (sd) {
	MATRIX(*sd, citingcat, xidx) +=
	  (xk-oldm)*(xk-MATRIX(*agd, citingcat, xidx));
      }
    }

    VECTOR(edges)[citingcat]+=igraph_vector_size(&neis);

    /* only the new node needs to be added */
    VECTOR(ntd)[citingcat]++;
    if (VECTOR(ntd)[citingcat]==1) {
      for (j=0;j<nocats;j++) {
	MATRIX(ch, j, citingcat)=VECTOR(edges)[j];
      }
    }
    
  }

  /* done, updates */
  for (i=0; i<nocats; i++) {
    for (j=0; j<nocats; j++) {
      igraph_real_t oldmean;
      if (VECTOR(ntd)[j] != 0) {
	MATRIX(normfact, i, j) += (VECTOR(edges)[i]-MATRIX(ch, i,  j)+1);
      }
      oldmean=MATRIX(*agd, i, j);
      MATRIX(*agd, i, j) *= MATRIX(notnull, i, j)/MATRIX(normfact,i, j);
      if (sd) {
	MATRIX(*sd, i, j) += oldmean * oldmean * MATRIX(notnull, i, j) *
	  (1-MATRIX(notnull, i, j)/MATRIX(normfact, i, j));
	if (MATRIX(normfact, i, j)>0) {
	  MATRIX(*sd, i, j)=
	    sqrt(MATRIX(*sd, i, j)/(MATRIX(normfact, i, j)-1));
	}
      }
    }
  }

  igraph_vector_destroy(&edges);
  igraph_vector_destroy(&neis);
  
  if (no) {
    igraph_matrix_destroy(no);
    *no=normfact;
  } else {
    igraph_matrix_destroy(&normfact);
  }
  
  igraph_matrix_destroy(&notnull);
  igraph_matrix_destroy(&ch);
  igraph_vector_destroy(&ntd);
  IGRAPH_FINALLY_CLEAN(6);
  return 0;
}
  
int igraph_measure_dynamics_citingcat_citedcat_st(const igraph_t *graph,
						  igraph_vector_t *res,
						  const igraph_matrix_t *agd,
						  const igraph_vector_t *cats,
						  igraph_integer_t pnocats) {
  long int no_of_nodes=igraph_vcount(graph);
  long int nocats=pnocats;
  
  igraph_matrix_t allst;
  long int j, node;

  IGRAPH_MATRIX_INIT_FINALLY(&allst, nocats, no_of_nodes);
  IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
  for (j=0; j<nocats; j++) {
    MATRIX(allst, j, 0)=MATRIX(*agd, j, (long int)VECTOR(*cats)[0]);
  }
  VECTOR(*res)[0]=MATRIX(allst, (long int) VECTOR(*cats)[0], 0);
  
  for (node=1; node<no_of_nodes; node++) {
    long int citingcat=VECTOR(*cats)[node];

    IGRAPH_ALLOW_INTERRUPTION();

    /* new node */
    for (j=0; j<nocats; j++) {
      MATRIX(allst, j, node)=MATRIX(allst, j, node-1) + 
	MATRIX(*agd, j, citingcat);
    }
    VECTOR(*res)[node]=MATRIX(allst, citingcat, node);
  }
  
  igraph_matrix_destroy(&allst);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}
