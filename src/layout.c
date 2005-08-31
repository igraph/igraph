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
#include "random.h"

int igraph_layout_random(igraph_t *graph, vector_t *res) {
  /* TODO */
  return 0;
}

int igraph_layout_circle(igraph_t *graph, vector_t *res) {
  /* TODO */
  return 0;
}

int igraph_layout_fruchterman_reingold(igraph_t *graph, vector_t *res, 
				       integer_t niter, real_t coolexp,
				       integer_t frame, vector_t *initial,
				       real_t initemp) {
  /* TODO */
  return 0;
}

int igraph_layout_kamada_kawai(igraph_t *graph, matrix_t *res,
			       integer_t niter, real_t sigma, 
			       real_t initemp, real_t coolexp,
			       real_t kkconst) {

  real_t temp, candx, candy;
  real_t dpot, odis, ndis, osqd, nsqd;
  long int n,i,j,k;
  matrix_t elen;
  vector_t vids;

  /* Define various things */
  n=igraph_vcount(graph);

  /* Calculate elen, initial x & y */

  RNG_BEGIN();

  matrix_resize(res, n, 2);
  matrix_init(&elen, n, n);
  vector_init_seq(&vids, 0, n-1);
  igraph_shortest_paths(graph, &elen, &vids, 3);
  vector_destroy(&vids);
  for (i=0; i<n; i++) {
    MATRIX(elen, i, i) = sqrt(n);
    MATRIX(*res, i, 0) = RNG_NORMAL(0, n/4.0);
    MATRIX(*res, i, 1) = RNG_NORMAL(0, n/4.0);
  }
  
  /*Perform the annealing loop*/
  temp=initemp;
  for(i=0;i<niter;i++){
    /*Update each vertex*/
    for(j=0;j<n;j++){
      /*Draw the candidate via a gaussian perturbation*/
      candx=RNG_NORMAL(MATRIX(*res, j, 0),sigma*temp/initemp);
      candy=RNG_NORMAL(MATRIX(*res, j, 1),sigma*temp/initemp);
      /*Calculate the potential difference for the new position*/
      dpot=0.0;
      for(k=0;k<n;k++)  /*Potential differences for pairwise effects*/
        if(j!=k){
          odis=sqrt((MATRIX(*res, j, 0)-MATRIX(*res, k, 0))*
		    (MATRIX(*res, j, 0)-MATRIX(*res, k, 0))+
		    (MATRIX(*res, j, 1)-MATRIX(*res, k, 1))*
		    (MATRIX(*res, j, 1)-MATRIX(*res, k, 1)));
          ndis=sqrt((candx-MATRIX(*res, k, 0))*(candx-MATRIX(*res, k, 0))+
		    (candy-MATRIX(*res, k, 1))*(candy-MATRIX(*res, k, 1)));
          osqd=(odis-MATRIX(elen, j, k))*(odis-MATRIX(elen, j, k));
          nsqd=(ndis-MATRIX(elen, j, k))*(ndis-MATRIX(elen, j, k));
          dpot+=kkconst*(osqd-nsqd)/(MATRIX(elen, j, k)*MATRIX(elen, j, k));
        }
      /*Make a keep/reject decision*/
      if(log(RNG_UNIF(0.0,1.0))<dpot/temp){
        MATRIX(*res, j, 0)=candx;
        MATRIX(*res, j, 1)=candy;
      }
    }
    /*Cool the system*/
    temp*=coolexp;
  }

  RNG_END();
  matrix_destroy(&elen);

  return 0;
}

int igraph_layout_springs(igraph_t *graph, vector_t *res,
			  real_t mass, real_t equil, real_t k,
			  real_t repeqdis, real_t kfr, bool_t repulse) {
  /* TODO */
  return 0;
}
