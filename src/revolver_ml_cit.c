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

  igraph_progress("ML Revolver d", 0, NULL);
  
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

    igraph_progress("ML Revolver d", 100*(it+1)/niter, NULL);

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
			 
