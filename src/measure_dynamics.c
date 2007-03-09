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

#include <math.h>

int igraph_measure_dynamics_id(const igraph_t *graph, igraph_integer_t start_vertex,
			       igraph_matrix_t *ak, igraph_matrix_t *sd,
			       igraph_matrix_t *confint, igraph_matrix_t *no,
			       const igraph_vector_t *st, igraph_integer_t pmaxind,
			       igraph_real_t significance, igraph_bool_t lno) {
  
  long int maxind=pmaxind;
  long int no_of_nodes=igraph_vcount(graph);

  int *indegree;
  igraph_matrix_t normfact;
  igraph_vector_t ntk, ch, notnull;
  igraph_vector_t neis;
  
  long int node;
  long int i;
  long int edges=0;
  
  igraph_bool_t lsd=(significance != 0);

  igraph_vector_init(&neis, 0);
  indegree=Calloc(no_of_nodes, int);
  igraph_matrix_resize(ak, maxind+1, 1);
  igraph_matrix_null(ak);
  if (lsd) {
    igraph_matrix_resize(sd, maxind+1, 1);
    igraph_matrix_resize(confint, maxind+1, 1);
    igraph_matrix_null(ak);
    /* TODO: significance calculation */
  }
  igraph_vector_init(&ntk, maxind+1);
  igraph_vector_init(&ch, maxind+1);
  igraph_matrix_init(&normfact, maxind+1, 1);
  igraph_vector_init(&notnull, maxind+1);

  for (node=0; node<start_vertex; node++) {
    
    igraph_neighbors(graph, &neis, node, IGRAPH_OUT);
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=indegree[to];
      
      indegree[to] ++;
      VECTOR(ntk)[xidx] --;
      VECTOR(ntk)[xidx+1] ++;
    }
    
    VECTOR(ntk)[0]++;
  }
  
  if (start_vertex != 0) {
    for (i=0; i<maxind+1; i++) {
      VECTOR(ch)[i]=start_vertex;
    }
  }
  
  for (node=start_vertex; node<no_of_nodes; node++) {

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
	MATRIX(*confint, xidx, 0) += (xk-oldm)*(xk-MATRIX(*ak, xidx, 0));
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
      /* TODO: confidence interval estimation */
      MATRIX(*confint, i, 0) +=
	oldmean * oldmean * VECTOR(notnull)[i] *
	(1-VECTOR(notnull)[i]/MATRIX(normfact, i, 0));
      if (MATRIX(normfact, i, 0) > 0) {
	MATRIX(*confint, i, 0) =
	  sqrt(MATRIX(*confint, i, 0)/(MATRIX(normfact, i, 0)-1));
	MATRIX(*sd, i, 0) = MATRIX(*confint,i,0);
      }
    }
  }
  
  if (!lno) {
    igraph_matrix_destroy(&normfact);
  } else {
    igraph_matrix_destroy(no);
    *no=normfact;
  }
  
  Free(indegree);
  igraph_vector_destroy(&ntk);
  igraph_vector_destroy(&ch);
  igraph_vector_destroy(&notnull);
  igraph_vector_destroy(&neis);
  
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
  
  indegree=Calloc(no_of_nodes, int);
  
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
  Free(indegree);

  return 0;
}

int igraph_measure_dynamics_idwindow(const igraph_t *graph, 
				     igraph_integer_t start_vertex,
				     igraph_matrix_t *ak, 
				     igraph_matrix_t *sd,
				     igraph_matrix_t *confint,
				     igraph_matrix_t *no,
				     const igraph_vector_t *st,
				     igraph_integer_t pmaxind,
				     igraph_real_t significance,
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
  
  igraph_bool_t lsd=(significance != 0);
  igraph_bool_t lno=(no != 0);
  
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
    igraph_matrix_resize(confint, maxind+1, 1);
  } 
  
  for (node=0; node<no_of_nodes; node++) {

    IGRAPH_ALLOW_INTERRUPTION();

    /* estimate Ak */
    igraph_neighbors(graph, &neis, node, IGRAPH_OUT);
    if (node>=start_vertex) {
      for (i=0; i<igraph_vector_size(&neis); i++) {
	long int to=VECTOR(neis)[i];
	long int xidx=VECTOR(indegree)[to];
	
	double xk=VECTOR(*st)[node]/VECTOR(ntk)[xidx];
	double oldm=MATRIX(*ak, xidx, 0);
	VECTOR(notnull)[xidx]+=1;
	MATRIX(*ak, xidx, 0) += (xk-oldm)/VECTOR(notnull)[xidx];
	if (lsd) {
	  MATRIX(*confint, xidx, 0) += (xk-oldm)*(xk-MATRIX(*ak, xidx, 0));
	}
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
      MATRIX(*confint, i, 0) +=
	oldmean * oldmean * VECTOR(notnull)[i] *
	(1-VECTOR(notnull)[i]/MATRIX(normfact, i, 0));
      if (MATRIX(normfact, i, 0) > 0) {
	MATRIX(*confint, i, 0) =
	  sqrt(MATRIX(*confint, i, 0)/(MATRIX(normfact, i, 0)-1));
	MATRIX(*sd, i, 0) = MATRIX(*confint,i,0);
      }
    }
  }
  
  igraph_dqueue_destroy(&history);
  igraph_vector_destroy(&notnull);
  if (!lno) {
    igraph_matrix_destroy(&normfact);
  } else {
    igraph_matrix_destroy(no);
    *no=normfact;
  }
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
    
  

/* alpha = 
    0.1,  0.05,   0.025,  0.01,   0.005,  0.001 */
const int igraph_i_tuppercrit_length=100;
const igraph_real_t igraph_i_tuppercrit[][6] = {
  { 3.078,  6.314, 12.706, 31.821, 63.657,318.313 }, /* dof: 1 */
  { 1.886,  2.920,  4.303,  6.965,  9.925, 22.327 },
  { 1.638,  2.353,  3.182,  4.541,  5.841, 10.215 },
  { 1.533,  2.132,  2.776,  3.747,  4.604,  7.173 },
  { 1.476,  2.015,  2.571,  3.365,  4.032,  5.893 },
  { 1.440,  1.943,  2.447,  3.143,  3.707,  5.208 },
  { 1.415,  1.895,  2.365,  2.998,  3.499,  4.782 },
  { 1.397,  1.860,  2.306,  2.896,  3.355,  4.499 },
  { 1.383,  1.833,  2.262,  2.821,  3.250,  4.296 },
  { 1.372,  1.812,  2.228,  2.764,  3.169,  4.143 },
  { 1.363,  1.796,  2.201,  2.718,  3.106,  4.024 },
  { 1.356,  1.782,  2.179,  2.681,  3.055,  3.929 },
  { 1.350,  1.771,  2.160,  2.650,  3.012,  3.852 },
  { 1.345,  1.761,  2.145,  2.624,  2.977,  3.787 },
  { 1.341,  1.753,  2.131,  2.602,  2.947,  3.733 },
  { 1.337,  1.746,  2.120,  2.583,  2.921,  3.686 },
  { 1.333,  1.740,  2.110,  2.567,  2.898,  3.646 },
  { 1.330,  1.734,  2.101,  2.552,  2.878,  3.610 },
  { 1.328,  1.729,  2.093,  2.539,  2.861,  3.579 },
  { 1.325,  1.725,  2.086,  2.528,  2.845,  3.552 },
  { 1.323,  1.721,  2.080,  2.518,  2.831,  3.527 },
  { 1.321,  1.717,  2.074,  2.508,  2.819,  3.505 },
  { 1.319,  1.714,  2.069,  2.500,  2.807,  3.485 },
  { 1.318,  1.711,  2.064,  2.492,  2.797,  3.467 },
  { 1.316,  1.708,  2.060,  2.485,  2.787,  3.450 },
  { 1.315,  1.706,  2.056,  2.479,  2.779,  3.435 },
  { 1.314,  1.703,  2.052,  2.473,  2.771,  3.421 },
  { 1.313,  1.701,  2.048,  2.467,  2.763,  3.408 },
  { 1.311,  1.699,  2.045,  2.462,  2.756,  3.396 },
  { 1.310,  1.697,  2.042,  2.457,  2.750,  3.385 },
  { 1.309,  1.696,  2.040,  2.453,  2.744,  3.375 },
  { 1.309,  1.694,  2.037,  2.449,  2.738,  3.365 },
  { 1.308,  1.692,  2.035,  2.445,  2.733,  3.356 },
  { 1.307,  1.691,  2.032,  2.441,  2.728,  3.348 },
  { 1.306,  1.690,  2.030,  2.438,  2.724,  3.340 },
  { 1.306,  1.688,  2.028,  2.434,  2.719,  3.333 },
  { 1.305,  1.687,  2.026,  2.431,  2.715,  3.326 },
  { 1.304,  1.686,  2.024,  2.429,  2.712,  3.319 },
  { 1.304,  1.685,  2.023,  2.426,  2.708,  3.313 },
  { 1.303,  1.684,  2.021,  2.423,  2.704,  3.307 },
  { 1.303,  1.683,  2.020,  2.421,  2.701,  3.301 },
  { 1.302,  1.682,  2.018,  2.418,  2.698,  3.296 },
  { 1.302,  1.681,  2.017,  2.416,  2.695,  3.291 },
  { 1.301,  1.680,  2.015,  2.414,  2.692,  3.286 },
  { 1.301,  1.679,  2.014,  2.412,  2.690,  3.281 },
  { 1.300,  1.679,  2.013,  2.410,  2.687,  3.277 },
  { 1.300,  1.678,  2.012,  2.408,  2.685,  3.273 },
  { 1.299,  1.677,  2.011,  2.407,  2.682,  3.269 },
  { 1.299,  1.677,  2.010,  2.405,  2.680,  3.265 },
  { 1.299,  1.676,  2.009,  2.403,  2.678,  3.261 },
  { 1.298,  1.675,  2.008,  2.402,  2.676,  3.258 },
  { 1.298,  1.675,  2.007,  2.400,  2.674,  3.255 },
  { 1.298,  1.674,  2.006,  2.399,  2.672,  3.251 },
  { 1.297,  1.674,  2.005,  2.397,  2.670,  3.248 },
  { 1.297,  1.673,  2.004,  2.396,  2.668,  3.245 },
  { 1.297,  1.673,  2.003,  2.395,  2.667,  3.242 },
  { 1.297,  1.672,  2.002,  2.394,  2.665,  3.239 },
  { 1.296,  1.672,  2.002,  2.392,  2.663,  3.237 },
  { 1.296,  1.671,  2.001,  2.391,  2.662,  3.234 },
  { 1.296,  1.671,  2.000,  2.390,  2.660,  3.232 },
  { 1.296,  1.670,  2.000,  2.389,  2.659,  3.229 },
  { 1.295,  1.670,  1.999,  2.388,  2.657,  3.227 },
  { 1.295,  1.669,  1.998,  2.387,  2.656,  3.225 },
  { 1.295,  1.669,  1.998,  2.386,  2.655,  3.223 },
  { 1.295,  1.669,  1.997,  2.385,  2.654,  3.220 },
  { 1.295,  1.668,  1.997,  2.384,  2.652,  3.218 },
  { 1.294,  1.668,  1.996,  2.383,  2.651,  3.216 },
  { 1.294,  1.668,  1.995,  2.382,  2.650,  3.214 },
  { 1.294,  1.667,  1.995,  2.382,  2.649,  3.213 },
  { 1.294,  1.667,  1.994,  2.381,  2.648,  3.211 },
  { 1.294,  1.667,  1.994,  2.380,  2.647,  3.209 },
  { 1.293,  1.666,  1.993,  2.379,  2.646,  3.207 },
  { 1.293,  1.666,  1.993,  2.379,  2.645,  3.206 },
  { 1.293,  1.666,  1.993,  2.378,  2.644,  3.204 },
  { 1.293,  1.665,  1.992,  2.377,  2.643,  3.202 },
  { 1.293,  1.665,  1.992,  2.376,  2.642,  3.201 },
  { 1.293,  1.665,  1.991,  2.376,  2.641,  3.199 },
  { 1.292,  1.665,  1.991,  2.375,  2.640,  3.198 },
  { 1.292,  1.664,  1.990,  2.374,  2.640,  3.197 },
  { 1.292,  1.664,  1.990,  2.374,  2.639,  3.195 },
  { 1.292,  1.664,  1.990,  2.373,  2.638,  3.194 },
  { 1.292,  1.664,  1.989,  2.373,  2.637,  3.193 },
  { 1.292,  1.663,  1.989,  2.372,  2.636,  3.191 },
  { 1.292,  1.663,  1.989,  2.372,  2.636,  3.190 },
  { 1.292,  1.663,  1.988,  2.371,  2.635,  3.189 },
  { 1.291,  1.663,  1.988,  2.370,  2.634,  3.188 },
  { 1.291,  1.663,  1.988,  2.370,  2.634,  3.187 },
  { 1.291,  1.662,  1.987,  2.369,  2.633,  3.185 },
  { 1.291,  1.662,  1.987,  2.369,  2.632,  3.184 },
  { 1.291,  1.662,  1.987,  2.368,  2.632,  3.183 },
  { 1.291,  1.662,  1.986,  2.368,  2.631,  3.182 },
  { 1.291,  1.662,  1.986,  2.368,  2.630,  3.181 },
  { 1.291,  1.661,  1.986,  2.367,  2.630,  3.180 },
  { 1.291,  1.661,  1.986,  2.367,  2.629,  3.179 },
  { 1.291,  1.661,  1.985,  2.366,  2.629,  3.178 },
  { 1.290,  1.661,  1.985,  2.366,  2.628,  3.177 },
  { 1.290,  1.661,  1.985,  2.365,  2.627,  3.176 },
  { 1.290,  1.661,  1.984,  2.365,  2.627,  3.175 },
  { 1.290,  1.660,  1.984,  2.365,  2.626,  3.175 },
  { 1.290,  1.660,  1.984,  2.364,  2.626,  3.174 }, /* dof: 100 */
  { 1.282,  1.645,  1.960,  2.326,  2.576,  3.090 }  /* dof: infinity */
};

int igraph_measure_dynamics_idage(const igraph_t *graph, igraph_integer_t start_vertex,
				  igraph_matrix_t *akl, 
				  igraph_matrix_t *sd, igraph_matrix_t *confint, 
				  igraph_matrix_t *no,
				  const igraph_vector_t *st, igraph_integer_t pagebins,
				  igraph_integer_t pmaxind, igraph_real_t significance,
				  igraph_bool_t lno) {

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

  igraph_bool_t lsd=(significance != 0);
  int signidx=0;

  binwidth = no_of_nodes/agebins+1;

  igraph_vector_init(&neis, 0);
  indegree=Calloc(no_of_nodes, int);
  igraph_matrix_resize(akl, maxind+1, agebins);
  igraph_matrix_null(akl);
  if (lsd) {
    igraph_matrix_resize(sd, maxind+1, agebins);
    igraph_matrix_resize(confint, maxind+1, agebins);
    igraph_matrix_null(sd);
    if (significance >= 0.998) {
      signidx=5;
    } else if (significance >= 0.99) {
      signidx=4;
    } else if (significance >= 0.98) {
      signidx=3;
    } else if (significance >= 0.95) {
      signidx=2;
    } else if (significance >= 0.9) {
      signidx=1;
    } else {
      signidx=0;
    }
  }
  igraph_matrix_init(&ntkl, maxind+1, agebins+1);
  igraph_matrix_init(&ch, maxind+1, agebins+1);
  igraph_matrix_init(&normfact, maxind+1, agebins);
  igraph_matrix_init(&notnull, maxind+1, agebins);
  
  for (node=0; node<start_vertex; node++) {

    MATRIX(ntkl, 0, (long int)(start_vertex-node)/binwidth) ++;

    igraph_neighbors(graph, &neis, node, IGRAPH_OUT);
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=indegree[to];
      long int yidx=(start_vertex-to)/binwidth;
      
      indegree[to] ++;
      MATRIX(ntkl, xidx, yidx)--;
      MATRIX(ntkl, xidx+1, yidx)++;
    }
  }

  if (start_vertex != 0) {
    for (i=0; i<maxind+1; i++) {
      for (j=0; j<agebins; j++) {
	MATRIX(ch, i, j) = start_vertex;
      }
    }
  }

  for (node=start_vertex; node<no_of_nodes; node++) {
    
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
	MATRIX(*confint, xidx, yidx) += (xk-oldm)*(xk-MATRIX(*akl, xidx, yidx));
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
    igraph_real_t tuppercrit;
    for (j=0; j<agebins; j++) {
      igraph_real_t oldmean;
      if (MATRIX(ntkl, i, j) != 0) {
	MATRIX(normfact, i, j) += (edges-MATRIX(ch, i, j)+1);
      }
      oldmean=MATRIX(*akl, i, j);
      MATRIX(*akl, i, j) *= MATRIX(notnull, i, j) / MATRIX(normfact, i, j);
      if (lsd) {
	MATRIX(*confint, i, j) +=
	  oldmean * oldmean * MATRIX(notnull, i, j) * 
	  (1-MATRIX(notnull,i,j)/MATRIX(normfact,i,j));
	if (MATRIX(normfact,i,j) > 0) {
	  MATRIX(*confint, i, j) =
	    sqrt(MATRIX(*confint, i, j)/(MATRIX(normfact,i,j)-1));
	  if (MATRIX(normfact,i,j) > igraph_i_tuppercrit_length) {
	    tuppercrit=
	      igraph_i_tuppercrit[igraph_i_tuppercrit_length][signidx];
	  } else {
	    tuppercrit=igraph_i_tuppercrit
	      [(long int)MATRIX(normfact,i,j)-1][signidx];
	  }	  
	  MATRIX(*sd, i, j) = MATRIX(*confint,i,j);
	  MATRIX(*confint, i, j) =
	    tuppercrit * MATRIX(*confint,i,j)/sqrt(MATRIX(normfact,i,j));
	}
      }
    }
  }
  
  if (!lno) {
    igraph_matrix_destroy(&normfact);
  } else {
    igraph_matrix_destroy(no);
    *no=normfact;
  }

  Free(indegree);
  igraph_matrix_destroy(&ntkl);
  igraph_matrix_destroy(&ch);
  igraph_matrix_destroy(&notnull);
  igraph_vector_destroy(&neis);

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
  
  indegree=Calloc(no_of_nodes, int);
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
  Free(indegree);
  
  return 0;
}

int igraph_measure_dynamics_idage_debug(const igraph_t *graph, igraph_matrix_t *akl,
					igraph_matrix_t *sd, igraph_matrix_t *confint, 
					igraph_matrix_t *no,
					const igraph_vector_t *st, igraph_integer_t pagebins,
					igraph_integer_t pmaxind, igraph_real_t significance,
					igraph_vector_t *estimates, 
					igraph_integer_t est_ind, igraph_integer_t est_age,
					igraph_bool_t lno) {
  
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

  igraph_bool_t lsd=(significance != 0);
  int signidx=0;

  binwidth = no_of_nodes/agebins+1;

  igraph_vector_init(&neis, 0);
  indegree=Calloc(no_of_nodes, int);
  igraph_matrix_resize(akl, maxind+1, agebins);
  igraph_matrix_null(akl);
  if (lsd) {
    igraph_matrix_resize(sd, maxind+1, agebins);
    igraph_matrix_null(sd);
    if (significance >= 0.998) {
      signidx=5;
    } else if (significance >= 0.99) {
      signidx=4;
    } else if (significance >= 0.98) {
      signidx=3;
    } else if (significance >= 0.95) {
      signidx=2;
    } else if (significance >= 0.9) {
      signidx=1;
    } else {
      signidx=0;
    }
  }

  igraph_vector_clear(estimates);

  igraph_matrix_init(&ntkl, maxind+1, agebins+1);
  igraph_matrix_init(&ch, maxind+1, agebins+1);
  igraph_matrix_init(&normfact, maxind+1, agebins);
  igraph_matrix_init(&notnull, maxind+1, agebins);

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
	MATRIX(*confint, xidx, yidx) += (xk-oldm)*(xk-MATRIX(*akl, xidx, yidx));
      }

      if (xidx==est_ind && yidx==est_age) {
	igraph_vector_push_back(estimates, xk);
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
    igraph_real_t tuppercrit;
    for (j=0; j<agebins; j++) {
      igraph_real_t oldmean;
      if (MATRIX(ntkl, i, j) != 0) {
	MATRIX(normfact, i, j) += (edges-MATRIX(ch, i, j)+1);
      }
      oldmean=MATRIX(*akl, i, j);
      MATRIX(*akl, i, j) *= MATRIX(notnull, i, j) / MATRIX(normfact, i, j);
      if (lsd) {
	MATRIX(*confint, i, j) +=
	  oldmean * oldmean * MATRIX(notnull, i, j) * 
	  (1-MATRIX(notnull,i,j)/MATRIX(normfact,i,j));
	if (MATRIX(normfact,i,j) > 0) {
	  MATRIX(*confint, i, j) =
	    sqrt(MATRIX(*confint, i, j)/(MATRIX(normfact,i,j)-1));
	  if (MATRIX(normfact,i,j) > igraph_i_tuppercrit_length) {
	    tuppercrit=
	      igraph_i_tuppercrit[igraph_i_tuppercrit_length][signidx];
	  } else {
	    tuppercrit=igraph_i_tuppercrit
	      [(long int)MATRIX(normfact,i,j)-1][signidx];
	  }
	  MATRIX(*sd, i, j) = MATRIX(*confint,i,j);
	  MATRIX(*confint, i, j) =
	    tuppercrit * MATRIX(*confint,i,j)/sqrt(MATRIX(normfact,i,j));
	}
      }
    }
  }

  for (i=0; i<igraph_vector_size(estimates); i++) {
    VECTOR(*estimates)[i] /= MATRIX(normfact, (long int)est_ind, 
				    (long int)est_age);
  }
/*   igraph_vector_push_back(estimates, MATRIX(normfact, est_ind, est_age)); */


  if (!lno) {
    igraph_matrix_destroy(&normfact);
  } else {
    igraph_matrix_destroy(no);
    *no=normfact;
  }
  
  Free(indegree);
  igraph_matrix_destroy(&ntkl);
  igraph_matrix_destroy(&ch);
  igraph_matrix_destroy(&notnull);
  igraph_vector_destroy(&neis);

  return 0;
}

int igraph_measure_dynamics_idwindowage(const igraph_t *graph, 
					igraph_integer_t start_vertex,
					igraph_matrix_t *akl, 
					igraph_matrix_t *sd, 
					igraph_matrix_t *confint, 
					igraph_matrix_t *no,
					const igraph_vector_t *st, 
					igraph_integer_t pagebins,
					igraph_integer_t pmaxind, 
					igraph_real_t significance,
					igraph_bool_t lno, 
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

  igraph_bool_t lsd=(significance != 0);
  int signidx=0;

  igraph_dqueue_t history;

  binwidth = no_of_nodes/agebins+1;

  igraph_vector_init(&neis, 0);
  indegree=Calloc(no_of_nodes, int);
  igraph_matrix_resize(akl, maxind+1, agebins);
  igraph_matrix_null(akl);
  if (lsd) {
    igraph_matrix_resize(sd, maxind+1, agebins);
    igraph_matrix_resize(confint, maxind+1, agebins);
    igraph_matrix_null(sd);
    if (significance >= 0.998) {
      signidx=5;
    } else if (significance >= 0.99) {
      signidx=4;
    } else if (significance >= 0.98) {
      signidx=3;
    } else if (significance >= 0.95) {
      signidx=2;
    } else if (significance >= 0.9) {
      signidx=1;
    } else {
      signidx=0;
    }
  }
  igraph_matrix_init(&ntkl, maxind+1, agebins+1);
  igraph_matrix_init(&ch, maxind+1, agebins+1);
  igraph_matrix_init(&normfact, maxind+1, agebins);
  igraph_matrix_init(&notnull, maxind+1, agebins);
  igraph_dqueue_init(&history, time_window);
  
  for (node=0; node<start_vertex; node++) {

    MATRIX(ntkl, 0, (long int)(start_vertex-node)/binwidth) ++;

    igraph_neighbors(graph, &neis, node, IGRAPH_OUT);
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int to=VECTOR(neis)[i];
      long int xidx=indegree[to];
      long int yidx=(start_vertex-to)/binwidth;
      
      indegree[to] ++;
      MATRIX(ntkl, xidx, yidx)--;
      MATRIX(ntkl, xidx+1, yidx)++;
    }
  }

  if (start_vertex != 0) {
    for (i=0; i<maxind+1; i++) {
      for (j=0; j<agebins; j++) {
	MATRIX(ch, i, j) = start_vertex;
      }
    }
  }

  for (i=0; i<=start_vertex; i++) {
    igraph_dqueue_push(&history, -1);
  }

  for (node=start_vertex; node<no_of_nodes; node++) {
    
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
	MATRIX(*confint, xidx, yidx) += (xk-oldm)*(xk-MATRIX(*akl, xidx, yidx));
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
    igraph_real_t tuppercrit;
    for (j=0; j<agebins; j++) {
      igraph_real_t oldmean;
      if (MATRIX(ntkl, i, j) != 0) {
	MATRIX(normfact, i, j) += (edges-MATRIX(ch, i, j)+1);
      }
      oldmean=MATRIX(*akl, i, j);
      MATRIX(*akl, i, j) *= MATRIX(notnull, i, j) / MATRIX(normfact, i, j);
      if (lsd) {
	MATRIX(*confint, i, j) +=
	  oldmean * oldmean * MATRIX(notnull, i, j) * 
	  (1-MATRIX(notnull,i,j)/MATRIX(normfact,i,j));
	if (MATRIX(normfact,i,j) > 0) {
	  MATRIX(*confint, i, j) =
	    sqrt(MATRIX(*confint, i, j)/(MATRIX(normfact,i,j)-1));
	  if (MATRIX(normfact,i,j) > igraph_i_tuppercrit_length) {
	    tuppercrit=
	      igraph_i_tuppercrit[igraph_i_tuppercrit_length][signidx];
	  } else {
	    tuppercrit=igraph_i_tuppercrit
	      [(long int)MATRIX(normfact,i,j)-1][signidx];
	  }	  
	  MATRIX(*sd, i, j) = MATRIX(*confint,i,j);
	  MATRIX(*confint, i, j) =
	    tuppercrit * MATRIX(*confint,i,j)/sqrt(MATRIX(normfact,i,j));
	}
      }
    }
  }
  
  if (!lno) {
    igraph_matrix_destroy(&normfact);
  } else {
    igraph_matrix_destroy(no);
    *no=normfact;
  }

  igraph_dqueue_destroy(&history);
  Free(indegree);
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
  
  indegree=Calloc(no_of_nodes, int);
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
  Free(indegree);
  
  return 0;
}

int igraph_measure_dynamics_citedcat_id_age(const igraph_t *graph,
					    igraph_integer_t start_vertex,
					    igraph_array3_t *adkl,
					    igraph_array3_t *sd,
					    igraph_array3_t *confint,
					    igraph_array3_t *no,
					    const igraph_vector_t *st,
					    const igraph_vector_t *cats,
					    igraph_integer_t pno_cats,
					    igraph_integer_t pagebins,
					    igraph_integer_t pmaxind,
					    igraph_real_t significance,
					    igraph_bool_t lno) {
  
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
  
  igraph_bool_t lsd=(significance != 0);
  int signidx=0;

  binwidth=no_of_nodes/agebins+1;
  
  igraph_vector_init(&neis, 0);
  indegree=Calloc(no_of_nodes, int);
  igraph_array3_resize(adkl, no_cats, maxind+1, agebins);
  igraph_array3_null(adkl);
  if (lsd) {
    igraph_array3_resize(sd, no_cats, maxind+1, agebins);
    igraph_array3_null(sd);
    igraph_array3_resize(confint, no_cats, maxind+1, agebins);
    igraph_array3_null(confint);
    if (significance >= 0.998) {
      signidx=5;
    } else if (significance >= 0.99) {
      signidx=4;
    } else if (significance >= 0.98) {
      signidx=3;
    } else if (significance >= 0.95) {
      signidx=2;
    } else if (significance >= 0.9) {
      signidx=1;
    } else {
      signidx=0;
    }    
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
	ARRAY3(*confint, cidx, xidx, yidx) += 
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
      igraph_real_t tuppercrit;
      for (j=0; j<agebins; j++) {
	igraph_real_t oldmean;
	if (ARRAY3(ntkl, k, i, j) !=0) {
	  ARRAY3(normfact, k, i, j) += (edges-ARRAY3(ch, k, i, j)+1);
	}
	oldmean=ARRAY3(*adkl, k, i, j);
	ARRAY3(*adkl, k, i, j) *=
	  ARRAY3(notnull, k, i, j)/ARRAY3(normfact, k, i, j);
	if (lsd) {
	  ARRAY3(*confint, k, i, j) += 
	    oldmean * oldmean * ARRAY3(notnull, k, i, j)*
	    (1-ARRAY3(notnull, k, i, j)/ARRAY3(normfact, k, i, j));
	  if (ARRAY3(normfact, k, i, j) > 0) {
	    ARRAY3(*confint, k, i, j)=
	      sqrt(ARRAY3(*confint, k, i, j)/(ARRAY3(normfact, k, i, j)-1));
	    if (ARRAY3(normfact, k, i, j) > igraph_i_tuppercrit_length) {
	      tuppercrit=
		igraph_i_tuppercrit[igraph_i_tuppercrit_length][signidx];
	    } else {
	      tuppercrit=igraph_i_tuppercrit
		[(long int)ARRAY3(normfact, k, i, j)-1][signidx];
	    }
	    ARRAY3(*sd, k, i, j)=ARRAY3(*confint, k, i, k);
	    ARRAY3(*confint, k, i, j) = 
	      tuppercrit*ARRAY3(*confint,k,i,j)/sqrt(ARRAY3(normfact,k,i,j));
	  }
	}
      }
    }
  }
  
  if (!lno) {
    igraph_array3_destroy(&normfact);
  } else {
    igraph_array3_destroy(no);
    *no=normfact;
  }
  
  Free(indegree);
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
  indegree=Calloc(no_of_nodes, int);
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
  Free(indegree);
  
  return 0;
}

int igraph_measure_dynamics_citingcat_id_age(const igraph_t *graph,
					     igraph_integer_t start_vertex,
					     igraph_array3_t *adkl,
					     igraph_array3_t *sd,
					     igraph_array3_t *confint,
					     igraph_array3_t *no,
					     const igraph_vector_t *st,
					     const igraph_vector_t *cats,
					     igraph_integer_t pno_cats,
					     igraph_integer_t pagebins,
					     igraph_integer_t pmaxind,
					     igraph_real_t significance,
					     igraph_bool_t lno) {
  
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
  
  igraph_bool_t lsd=(significance != 0);
  int signidx=0;
  
  binwidth=no_of_nodes/agebins+1;
  
  igraph_vector_init(&neis, 0);
  indegree=Calloc(no_of_nodes, int);
  igraph_vector_init(&edges, no_cats);

  igraph_array3_resize(adkl, no_cats, maxind+1, agebins);
  igraph_array3_null(adkl);
  if (lsd) {
    igraph_array3_resize(sd, no_cats, maxind+1, agebins);
    igraph_array3_null(sd);
    igraph_array3_resize(confint, no_cats, maxind+1, agebins);
    igraph_array3_null(confint);
    if (significance >= 0.998) {
      signidx=5;
    } else if (significance >= 0.99) {
      signidx=4;
    } else if (significance >= 0.98) {
      signidx=3;
    } else if (significance >= 0.95) {
      signidx=2;
    } else if (significance >= 0.9) {
      signidx=1;
    } else {
      signidx=0;
    }    
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
	ARRAY3(*confint, cidx, xidx, yidx) += 
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
	}
	
      }
    }
  }
  
  if (!lno) {
    igraph_array3_destroy(&normfact);
  } else {
    igraph_array3_destroy(no);
    *no=normfact;
  }
  
  Free(indegree);
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
  indegree=Calloc(no_of_nodes, int);
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
  Free(indegree);
  
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
				igraph_matrix_t *no,
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
  IGRAPH_CHECK(igraph_vector_order(ntime, &ntimeidx, events));
  IGRAPH_CHECK(igraph_vector_order(etime, &etimeidx, events));

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
      MATRIX(notnull, yidx, xidx) += 1;

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
	  MATRIX(*sd, i, j) = sqrt(MATRIX(*sd, i, j)/MATRIX(normfact, i, j));
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
  if (!no) { 
    igraph_matrix_destroy(&normfact);
  } else {
    igraph_matrix_destroy(no);
    *no=normfact;
  }
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
  IGRAPH_CHECK(igraph_vector_order(ntime, &ntimeidx, events));
  IGRAPH_CHECK(igraph_vector_order(etime, &etimeidx, events));

  for (timestep=0; timestep<events; timestep++) {

    IGRAPH_ALLOW_INTERRUPTION();
    
    /* add the new nodes, if any */
    while (nptr < no_of_nodes &&
	   VECTOR(*ntime)[ (long int) VECTOR(ntimeidx)[nptr] ] == timestep) {
      igraph_real_t akk_inc;
      akk_inc=0;
      for (i=0; i<=maxdeg; i++) {
	akk_inc += VECTOR(ntk)[i]*MATRIX(*akk, i, 0);
	akk_inc += VECTOR(ntk)[i]*MATRIX(*akk, i, 0);
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
      akk_inc += MATRIX(*akk, yidx, yidx);
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
    
