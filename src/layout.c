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
#include <math.h>

/**
 * \section about_layouts
 * 
 * <para>Layout generator functions (or at least most of them) try place the
 * vertices and edges of a graph on a 2D plane or in 3D space in a way
 * which visually pleases the human eye.</para>
 *
 * <para>They take a graph object and a number of parameters as arguments
 * and return a <type>matrix_t</type>, in which each row gives the
 * coordinates of a vertex.</para>
 */

/**
 * \ingroup layout
 * \function igraph_layout_random
 * \brief Places the vertices uniform randomly on a plane.
 * 
 * \param graph Pointer to an initialized graph object.
 * \param res Pointer to an initialized graph object. This will
 *        contain the result and will be resized in needed.
 * \return Error code. The current implementation always returns with
 * success. 
 * 
 * Time complexity: O(|V|), the
 * number of vertices. 
 */

int igraph_layout_random(const igraph_t *graph, matrix_t *res) {

  long int no_of_nodes=igraph_vcount(graph);
  long int i;

  IGRAPH_CHECK(matrix_resize(res, no_of_nodes, 2));

  RNG_BEGIN();

  for (i=0; i<no_of_nodes; i++) {
    MATRIX(*res, i, 0)=RNG_UNIF(-1, 1);
    MATRIX(*res, i, 1)=RNG_UNIF(-1, 1);
  }

  RNG_END();
  
  return 0;
}

/**
 * \ingroup layout
 * \function igraph_layout_circle
 * \brief Places the vertices uniformly on a circle, in the order of
 * vertex ids.
 * 
 * \param graph Pointer to an initialized graph object.
 * \param res Pointer to an initialized graph object. This will
 *        contain the result and will be resized in needed.
 * \return Error code.
 * 
 * Time complexity: O(|V|), the
 * number of vertices. 
 */

int igraph_layout_circle(const igraph_t *graph, matrix_t *res) {

  long int no_of_nodes=igraph_vcount(graph);
  long int i;

  IGRAPH_CHECK(matrix_resize(res, no_of_nodes, 2));  

  for (i=0; i<no_of_nodes; i++) {
    real_t phi=2*M_PI/no_of_nodes*i;
    MATRIX(*res, i, 0)=cos(phi);
    MATRIX(*res, i, 1)=sin(phi);
  }
  
  return 0;
}

/**
 * \ingroup layout
 * \function igraph_layout_fruchterman_reingold
 * \brief Places the vertices on a plane according to the
 * Fruchterman-Reingold algorithm.
 *
 * This is a force-directed layout, see Fruchterman, T.M.J. and
 * Reingold, E.M.: Graph Drawing by Force-directed Placement.
 * Software -- Practice and Experience, 21/11, 1129--1164,
 * 1991. 
 * This function was ported from the SNA R package.
 * \param graph Pointer to an initialized graph object.
 * \param res Pointer to an initialized matrix object. This will
 *        contain the result and will be resized in needed.
 * \param niter The number of iterations to do.
 * \param maxdelta The maximum distance to move a vertex in an
 *        iteration.
 * \param area The area parameter of the algorithm.
 * \param coolexp The cooling exponent of the simulated annealing.
 * \param repulserad Determines the radius at which
 *        vertex-vertex repulsion cancels out attraction of
 *        adjacent vertices.
 * \param use_seed Logical, if true the supplied values in the
 *        <parameter>res</parameter> argument are used as an initial layout, if
 *        false a random initial layout is used.
 * \return Error code.
 * 
 * Time complexity: O(|V|^2) in each
 * iteration, |V| is the number of
 * vertices in the graph. 
 */

int igraph_layout_fruchterman_reingold(const igraph_t *graph, matrix_t *res,
				       integer_t niter, real_t maxdelta,
				       real_t area, real_t coolexp, 
				       real_t repulserad, bool_t use_seed) {
  real_t frk,t,ded,xd,yd;
  real_t rf,af;
  long int i,j,k;

  long int no_of_nodes=igraph_vcount(graph);
  matrix_t dxdy=MATRIX_NULL;
  igraph_es_t edgeit;
  
  IGRAPH_CHECK(matrix_resize(res, no_of_nodes, 2));
  if (!use_seed) {
    IGRAPH_CHECK(igraph_layout_random(graph, res));
  }
  MATRIX_INIT_FINALLY(&dxdy, no_of_nodes, 2);
  
  frk=sqrt(area/no_of_nodes);

  for(i=niter;i>0;i--) {
    /* Set the temperature (maximum move/iteration) */
    t=maxdelta*pow(i/(double)niter,coolexp);
    /* Clear the deltas */
    matrix_null(&dxdy);
    /* Increment deltas for each undirected pair */
    for(j=0;j<no_of_nodes;j++) {
      for(k=j+1;k<no_of_nodes;k++){
        /* Obtain difference vector */
        xd=MATRIX(*res, j, 0)-MATRIX(*res, k, 0);
        yd=MATRIX(*res, j, 1)-MATRIX(*res, k, 1);
        ded=sqrt(xd*xd+yd*yd);  /* Get dyadic euclidean distance */
        xd/=ded;                /* Rescale differences to length 1 */
        yd/=ded;
        /* Calculate repulsive "force" */
        rf=frk*frk*(1.0/ded-ded*ded/repulserad);
        MATRIX(dxdy, j, 0)+=xd*rf; /* Add to the position change vector */
        MATRIX(dxdy, k, 0)-=xd*rf;
        MATRIX(dxdy, j, 1)+=yd*rf;
        MATRIX(dxdy, k, 1)-=yd*rf;
      }
    }
    /* Calculate the attractive "force" */
    IGRAPH_CHECK(igraph_es_all(graph, &edgeit));
    while (!igraph_es_end(graph, &edgeit)) {
      j=igraph_es_from(graph, &edgeit);
      k=igraph_es_to(graph, &edgeit);
      xd=MATRIX(*res, j, 0)-MATRIX(*res, k, 0);
      yd=MATRIX(*res, j, 1)-MATRIX(*res, k, 1);
      ded=sqrt(xd*xd+yd*yd);  /* Get dyadic euclidean distance */
      xd/=ded;                /* Rescale differences to length 1 */
      yd/=ded;
      af=ded*ded/frk;
      MATRIX(dxdy, j, 0)-=xd*af; /* Add to the position change vector */
      MATRIX(dxdy, k, 0)+=xd*af;
      MATRIX(dxdy, j, 1)-=yd*af;
      MATRIX(dxdy, k, 1)+=yd*af;
      igraph_es_next(graph, &edgeit);
    }
    igraph_es_destroy(&edgeit);
    
    /* Dampen motion, if needed, and move the points */   
    for(j=0;j<no_of_nodes;j++){
      ded=sqrt(MATRIX(dxdy, j, 0)*MATRIX(dxdy, j, 0)+
	       MATRIX(dxdy, j, 1)*MATRIX(dxdy, j, 1));    
      if(ded>t){		/* Dampen to t */
        ded=t/ded;
        MATRIX(dxdy, j, 0)*=ded;
        MATRIX(dxdy, j, 1)*=ded;
      }
      MATRIX(*res, j, 0)+=MATRIX(dxdy, j, 0); /* Update positions */
      MATRIX(*res, j, 1)+=MATRIX(dxdy, j, 1);
    }
  }
  
  matrix_destroy(&dxdy);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

/**
 * \ingroup layout
 * \function igraph_layout_kamada_kawai
 * \brief Places the vertices on a plane according the Kamada-Kawai
 * algorithm. 
 *
 * This is a force directed layout, see  Kamada, T. and Kawai, S.: An
 * Algorithm for Drawing General Undirected Graphs. Information
 * Processing Letters, 31/1, 7--15, 1989.
 * This function was ported from the SNA R package.
 * \param graph A graph object.
 * \param res Pointer to an initialized matrix object. This will
 *        contain the result and will be resized if needed.
 * \param niter The number of iterations to perform.
 * \param sigma Sets the base standard deviation of position
 *        change proposals. 
 * \param initemp Sets the initial temperature for the annealing.
 * \param coolexp The cooling exponent of the annealing.
 * \param kkconst The Kamada-Kawai vertex attraction constant.
 * \return Error code.
 * 
 * Time complexity: O(|V|^2) for each
 * iteration, |V| is the number of
 * vertices in the graph. 
 */

int igraph_layout_kamada_kawai(const igraph_t *graph, matrix_t *res,
			       integer_t niter, real_t sigma, 
			       real_t initemp, real_t coolexp,
			       real_t kkconst) {

  real_t temp, candx, candy;
  real_t dpot, odis, ndis, osqd, nsqd;
  long int n,i,j,k;
  matrix_t elen;

  /* Define various things */
  n=igraph_vcount(graph);

  /* Calculate elen, initial x & y */

  RNG_BEGIN();

  IGRAPH_CHECK(matrix_resize(res, n, 2));
  MATRIX_INIT_FINALLY(&elen, n, n);
  IGRAPH_CHECK(igraph_shortest_paths(graph, &elen, IGRAPH_VS_ALL(graph), 3));
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
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

int igraph_layout_springs(const igraph_t *graph, matrix_t *res,
			  real_t mass, real_t equil, real_t k,
			  real_t repeqdis, real_t kfr, bool_t repulse) {
  /* TODO */
  return 0;
}
