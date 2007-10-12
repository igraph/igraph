/* -*- mode: C -*-  */
/* 
   IGraph R package.
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

#include "igraph.h"

#include <math.h>

typedef struct igraph_i_revolver_ml_D_data_t {
  const igraph_t *graph;
  long int no_of_nodes;
  igraph_revolver_ml_fdf_t *fdf;
  igraph_vector_t *A;
  igraph_vector_t *dA;
  igraph_integer_t maxdegree;
  igraph_vector_long_t *ptk;
  igraph_vector_t *neis;
  igraph_vector_long_t *degree;
} igraph_i_revolver_ml_D_data_t;  
    
igraph_real_t igraph_i_revolver_ml_D(igraph_real_t arg,
				      void *info) {
  
  igraph_i_revolver_ml_D_data_t *data=info;
  igraph_real_t sum=0;
  long int t;
  igraph_real_t S1=0, S2=0;
  
  fprintf(stderr, "Evaluating %f: ", (double)arg);
  
  /* Init */
  igraph_vector_long_null(data->ptk);
  igraph_vector_long_null(data->degree);

  /* Calculate all possible A and dA values */
  for (t=0; t<=data->maxdegree; t++) {
    igraph_real_t deg=t;
    data->fdf(&arg, &deg, VECTOR(*(data->A))+t, VECTOR(*(data->dA))+t);
  }

  for (t=0; t<data->no_of_nodes; t++) {
    long int n, nneis;
    IGRAPH_CHECK(igraph_neighbors(data->graph, data->neis, t, IGRAPH_OUT));
    nneis=igraph_vector_size(data->neis);

    /* Update sum */
    for (n=0; n<nneis; n++) {
      long int to=VECTOR(*(data->neis))[n];
      long int x=VECTOR(*(data->degree))[to];
      
      sum += VECTOR(*(data->dA))[x] / VECTOR(*(data->A))[x];
      sum -= S1/S2;
    }

    /* Update ptk */
    for (n=0; n<nneis; n++) {
      long int to=VECTOR(*(data->neis))[n];
      long int x=VECTOR(*(data->degree))[to];
      
      VECTOR(*(data->degree))[to] += 1;
      VECTOR(*(data->ptk))[x+1] += 1;
      VECTOR(*(data->ptk))[x] -= 1;
      S1 += VECTOR(*(data->dA))[x+1];
      S1 -= VECTOR(*(data->dA))[x];
      S2 += VECTOR(*(data->A))[x+1];
      S2 -= VECTOR(*(data->A))[x];
    }
    
    VECTOR(*(data->ptk))[0] += 1;
    S1 += VECTOR(*(data->dA))[0];
    S2 += VECTOR(*(data->A))[0];
  }

  fprintf(stderr, "%f\n", (double) sum);
  
  return sum;
}
			   

int igraph_revolver_ml_D(const igraph_t *graph,
			 igraph_real_t *res,
			 igraph_real_t left,
			 igraph_real_t right,
			 igraph_real_t delta,
			 int maxit,
			 igraph_revolver_ml_fdf_t *fdf) {
  
  igraph_i_revolver_ml_D_data_t info;
  long int no_of_nodes=igraph_vcount(graph);
  igraph_vector_t A, dA, neis;
  igraph_vector_long_t ptk, degree;
  igraph_integer_t maxdegree;
  int ret;

  IGRAPH_CHECK(igraph_maxdegree(graph, &maxdegree, igraph_vss_all(),
				IGRAPH_IN, IGRAPH_LOOPS));
  
  IGRAPH_CHECK(igraph_vector_long_init(&ptk, maxdegree+1));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &ptk);
  IGRAPH_CHECK(igraph_vector_long_init(&degree, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &degree);
  
  IGRAPH_VECTOR_INIT_FINALLY(&A, maxdegree+1);
  IGRAPH_VECTOR_INIT_FINALLY(&dA, maxdegree+1);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  
  info.graph=        graph;
  info.no_of_nodes=  no_of_nodes;
  info.fdf=          fdf;
  info.A=           &A;
  info.dA=          &dA;
  info.maxdegree=    maxdegree;
  info.ptk=         &ptk;
  info.neis=        &neis;
  info.degree=      &degree;
  
  /* Ok, call zero finding */
  ret=igraph_zeroin(left, right, igraph_i_revolver_ml_D, &info,
		    &delta, &maxit, res);

  igraph_vector_destroy(&neis);
  igraph_vector_destroy(&dA);
  igraph_vector_destroy(&A);
  igraph_vector_long_destroy(&degree);
  igraph_vector_long_destroy(&ptk);
  IGRAPH_FINALLY_CLEAN(5);

  if (ret != 0) {
    IGRAPH_WARNING("Iteration did not converge!");
  }  
  
  return 0;
}

int igraph_i_revolver_ml_D_alpha(const igraph_real_t *param,
				 const igraph_real_t *arg,
				 igraph_real_t *fres,
				 igraph_real_t *dfres) {

  if (*arg!=0) {
    igraph_real_t p=pow(*arg, *param);
    *fres = p+1.0;
    *dfres = p * log(*arg);
  } else {
    *fres = 1.0;
    *dfres= 0.0;
  }

  return 0;
}

int igraph_revolver_ml_D_alpha(const igraph_t *graph,
			       igraph_real_t *res,
			       igraph_real_t left,
			       igraph_real_t right,
			       igraph_real_t delta,
			       int maxit) {

  return igraph_revolver_ml_D(graph, res, left, right, delta, maxit, 
			      igraph_i_revolver_ml_D_alpha);
}

/* TODO: delta, logprob, lognull */

int igraph_revolver_ml_d(const igraph_t *graph,
			 igraph_integer_t niter,
			 igraph_vector_t *kernel,
			 igraph_vector_t *cites) {
  
  long int no_of_nodes=igraph_vcount(graph);
  igraph_integer_t imaxdegree;
  long int maxdegree, actmaxdegree;
  long int it, t, i;
  igraph_vector_long_t ptk;
  igraph_vector_t *mycites, vmycites;
  igraph_vector_t neis;
  igraph_vector_long_t degree;
  igraph_real_t S=0;

  igraph_vector_t vmykernel;
  igraph_vector_t *kernels[]={ kernel, &vmykernel };
  long int actkernel=0;
  igraph_vector_t *fromkernel=kernels[actkernel], 
    *tokernel=kernels[1-actkernel];
  
  IGRAPH_CHECK(igraph_maxdegree(graph, &imaxdegree, igraph_vss_all(),
				IGRAPH_IN, IGRAPH_LOOPS));
  maxdegree=imaxdegree;
  
  IGRAPH_CHECK(igraph_vector_long_init(&ptk, maxdegree+1));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &ptk);
  
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_CHECK(igraph_vector_long_init(&degree, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &degree);
  IGRAPH_VECTOR_INIT_FINALLY(&vmykernel, maxdegree+1);
  
  if (cites) {
    IGRAPH_CHECK(igraph_vector_resize(cites, maxdegree+1));
    igraph_vector_null(cites);
    mycites=cites;
  } else {
    IGRAPH_VECTOR_INIT_FINALLY(&vmycites, maxdegree+1);
    mycites=&vmycites;
  }

  IGRAPH_CHECK(igraph_vector_resize(kernel, maxdegree+1));
  igraph_vector_fill(kernel, 1);

  IGRAPH_PROGRESS("ML Revolver d", 0, NULL);
  
  for (it=0; it<niter; it++) {
    
    igraph_vector_null(tokernel);
    igraph_vector_long_null(&ptk);
    igraph_vector_long_null(&degree);
    S=0.0;
    actmaxdegree=0;

    for (t=0; t<no_of_nodes; t++) {
      long int n, nneis;
      IGRAPH_CHECK(igraph_neighbors(graph, &neis, t, IGRAPH_OUT));
      nneis=igraph_vector_size(&neis);

      /* Calculate some terms of the sum for the non-zero classes */
      if (S != 0) {
	for (i=0; i<=actmaxdegree; i++) {
	  VECTOR(*tokernel)[i] += nneis * VECTOR(ptk)[i] / S;
	}
      }
      
      /* Update ptk for the next time step */
      for (n=0; n<nneis; n++) {
	long int to=VECTOR(neis)[n];
	long int x=VECTOR(degree)[to];

	VECTOR(degree)[to] += 1;
	if (x==actmaxdegree) { actmaxdegree++; }

	VECTOR(ptk)[x+1] += 1;
	VECTOR(ptk)[x] -= 1;
	S += VECTOR(*fromkernel)[x+1];
	S -= VECTOR(*fromkernel)[x];
      
	if (it==0) {
	  VECTOR(*mycites)[x] += 1;
	}
      }
      VECTOR(ptk)[0] += 1;
      S += VECTOR(*fromkernel)[0];

    } /* t<no_of_nodes  */
    
    /* final step, Mk/sum */
    for (i=0; i<maxdegree; i++) {
      VECTOR(*tokernel)[i] = VECTOR(*mycites)[i] / VECTOR(*tokernel)[i];      
    }
    
    /* Switch kernels */
    actkernel=1-actkernel;
    fromkernel=kernels[actkernel];
    tokernel=kernels[1-actkernel];

    IGRAPH_PROGRESS("ML Revolver d", 100*(it+1)/niter, NULL);

  } /* it<niter */

  /* switch kernels if needed */
  if (fromkernel != kernel) {
    igraph_vector_clear(kernel);
    igraph_vector_append(kernel, fromkernel);
  }
  VECTOR(*kernel)[maxdegree]=IGRAPH_NAN;

  if (!cites) {
    igraph_vector_destroy(&vmycites);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  igraph_vector_destroy(&vmykernel);
  igraph_vector_long_destroy(&degree);
  igraph_vector_destroy(&neis);
  igraph_vector_long_destroy(&ptk);
  IGRAPH_FINALLY_CLEAN(4);
  return 0;
}
			 
