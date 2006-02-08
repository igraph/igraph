/* -*- mode: C -*-  */
/* 
   IGraph R package.
   Copyright (C) 2003, 2004, 2005, 2006  Gabor Csardi <csardi@rmki.kfki.hu>
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
 * <para>Layout generator functions (or at least most of them) try to place the
 * vertices and edges of a graph on a 2D plane or in 3D space in a way
 * which visually pleases the human eye.</para>
 *
 * <para>They take a graph object and a number of parameters as arguments
 * and return an \type igraph_matrix_t, in which each row gives the
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

int igraph_layout_random(const igraph_t *graph, igraph_matrix_t *res) {

  long int no_of_nodes=igraph_vcount(graph);
  long int i;

  IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 2));

  RNG_BEGIN();

  for (i=0; i<no_of_nodes; i++) {
    MATRIX(*res, i, 0)=RNG_UNIF(-1, 1);
    MATRIX(*res, i, 1)=RNG_UNIF(-1, 1);
  }

  RNG_END();
  
  return 0;
}

/**
 * \function igraph_layout_random_3d
 * \brief Random layout in 3D
 * 
 * \param graph The graph to place.
 * \param res Pointer to an initialized matrix object. It will be
 * resized to hold the result.
 * \return Error code. The current implementation always returns with
 * success. 
 *
 * Added in version 0.2.</para><para>
 * 
 * Time complexity: O(|V|), the number of vertices.
 */

int igraph_layout_random_3d(const igraph_t *graph, igraph_matrix_t *res) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int i;
  
  IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 3));
  
  RNG_BEGIN();
  
  for (i=0; i<no_of_nodes; i++) {
    MATRIX(*res, i, 0)=RNG_UNIF(-1, 1);
    MATRIX(*res, i, 1)=RNG_UNIF(-1, 1);
    MATRIX(*res, i, 2)=RNG_UNIF(-1, 1);
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

int igraph_layout_circle(const igraph_t *graph, igraph_matrix_t *res) {

  long int no_of_nodes=igraph_vcount(graph);
  long int i;

  IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 2));  

  for (i=0; i<no_of_nodes; i++) {
    real_t phi=2*M_PI/no_of_nodes*i;
    MATRIX(*res, i, 0)=cos(phi);
    MATRIX(*res, i, 1)=sin(phi);
  }
  
  return 0;
}

/**
 * \function igraph_layout_sphere
 * \brief Places vertices (more or less) uniformly on a sphere.
 * 
 * The algorithm was described in the following paper:
 * Distributing many points on a sphere by E.B. Saff and
 * A.B.J. Kuijlaars, \emb Mathematical Intelligencer \eme 19.1 (1997)
 * 5--11.  
 * 
 * \param graph Pointer to an initialized graph object.
 * \param res Pointer to an initialized matrix object, the will be
 * stored here. It will be resized.
 * \return Error code. The current implementation always returns with
 * success. 
 * 
 * Added in version 0.2.</para><para>
 *
 * Time complexity: O(|V|), the number of vertices in the graph.
 */

int igraph_layout_sphere(const igraph_t *graph, igraph_matrix_t *res) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int i;
  real_t h;
  
  IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 3));

  if (no_of_nodes != 0) {
    MATRIX(*res, 0, 0)=M_PI;
    MATRIX(*res, 0, 1)=0;
  }
  for (i=1; i<no_of_nodes-1; i++) {
    h = -1 + 2*i/(double)(no_of_nodes-1);
    MATRIX(*res, i, 0) = acos(h);
    MATRIX(*res, i, 1) = fmod((MATRIX(*res, i-1, 1) +
			       3.6/sqrt(no_of_nodes*(1-h*h))), 2*M_PI);
  }
  if (no_of_nodes >=2) {
    MATRIX(*res, no_of_nodes-1, 0)=0;
    MATRIX(*res, no_of_nodes-1, 1)=0;
  }

  for (i=0; i<no_of_nodes; i++) {
    real_t x=cos(MATRIX(*res, i, 1))*sin(MATRIX(*res, i, 0));
    real_t y=sin(MATRIX(*res, i, 1))*sin(MATRIX(*res, i, 0));
    real_t z=cos(MATRIX(*res, i, 0));
    MATRIX(*res, i, 0)=x;
    MATRIX(*res, i, 1)=y;
    MATRIX(*res, i, 2)=z;
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
 *        \p res argument are used as an initial layout, if
 *        false a random initial layout is used.
 * \return Error code.
 * 
 * Time complexity: O(|V|^2) in each
 * iteration, |V| is the number of
 * vertices in the graph. 
 */

int igraph_layout_fruchterman_reingold(const igraph_t *graph, igraph_matrix_t *res,
				       integer_t niter, real_t maxdelta,
				       real_t area, real_t coolexp, 
				       real_t repulserad, bool_t use_seed) {
  real_t frk,t,ded,xd,yd;
  real_t rf,af;
  long int i,j,k;

  long int no_of_nodes=igraph_vcount(graph);
  igraph_matrix_t dxdy=IGRAPH_MATRIX_NULL;
  igraph_es_t edgeit;
  
  IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 2));
  if (!use_seed) {
    IGRAPH_CHECK(igraph_layout_random(graph, res));
  }
  IGRAPH_MATRIX_INIT_FINALLY(&dxdy, no_of_nodes, 2);
  
  frk=sqrt(area/no_of_nodes);

  for(i=niter;i>0;i--) {
    /* Set the temperature (maximum move/iteration) */
    t=maxdelta*pow(i/(double)niter,coolexp);
    /* Clear the deltas */
    igraph_matrix_null(&dxdy);
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
      if (ded != 0) {
	xd/=ded;                /* Rescale differences to length 1 */
	yd/=ded;
      }
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
  
  igraph_matrix_destroy(&dxdy);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}

/**
 * \function igraph_layout_fruchterman_reingold_3d
 * \brief This is the 3D version of the force based
 * Fruchterman-Reingold layout (see \ref
 * igraph_layout_fruchterman_reingold for the 2D version
 *
 * This function was ported from the SNA R package.
 * \param graph Pointer to an initialized graph object.
 * \param res Pointer to an initialized matrix object. This will
 *        contain the result and will be resized in needed.
 * \param niter The number of iterations to do.
 * \param maxdelta The maximum distance to move a vertex in an
 *        iteration.
 * \param volume The volume parameter of the algorithm.
 * \param coolexp The cooling exponent of the simulated annealing.
 * \param repulserad Determines the radius at which
 *        vertex-vertex repulsion cancels out attraction of
 *        adjacent vertices.
 * \param use_seed Logical, if true the supplied values in the
 *        \p res argument are used as an initial layout, if
 *        false a random initial layout is used.
 * \return Error code.
 *
 * Added in version 0.2.</para><para>
 *
 * Time complexity: O(|V|^2) in each
 * iteration, |V| is the number of
 * vertices in the graph. 
 * 
 */

int igraph_layout_fruchterman_reingold_3d(const igraph_t *graph, 
					  igraph_matrix_t *res,
					  integer_t niter, real_t maxdelta,
					  real_t volume, real_t coolexp,
					  real_t repulserad,
					  bool_t use_seed) {
  
  real_t frk, t, ded, xd, yd, zd;
  igraph_matrix_t dxdydz;
  real_t rf, af;
  long int i, j, k;
  
  long int no_of_nodes=igraph_vcount(graph);
  igraph_es_t edgeit;
  
  IGRAPH_CHECK(igraph_matrix_init(&dxdydz, no_of_nodes, 3));
  IGRAPH_FINALLY(igraph_matrix_destroy, &dxdydz);
  
  IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 3));
  if (!use_seed) {
    IGRAPH_CHECK(igraph_layout_random_3d(graph, res));
  }

  frk=pow(volume/(double)no_of_nodes,1.0/3.0); /*Define the F-R constant*/

  /*Run the annealing loop*/
  for(i=niter;i>=0;i--){
    /*Set the temperature (maximum move/iteration)*/
    t=maxdelta*pow(i/(double)niter,coolexp);
    /*Clear the deltas*/
    igraph_matrix_null(&dxdydz);
    /*Increment deltas for each undirected pair*/
    for(j=0;j<no_of_nodes;j++) {
      for(k=j+1;k<no_of_nodes;k++){
        /*Obtain difference vector*/
        xd=MATRIX(*res, j, 0)-MATRIX(*res, k, 0);
        yd=MATRIX(*res, j, 1)-MATRIX(*res, k, 1);
        zd=MATRIX(*res, j, 2)-MATRIX(*res, k, 2);
        ded=sqrt(xd*xd+yd*yd+zd*zd);  /*Get dyadic euclidean distance*/
        if (ded != 0) {
	  xd/=ded;                      /*Rescale differences to length 1*/
	  yd/=ded;
	  zd/=ded;
	}
        /*Calculate repulsive "force"*/
        rf=frk*frk*(1.0/ded-ded*ded/repulserad);
        MATRIX(dxdydz, j, 0)+=xd*rf;     /*Add to the position change vector*/
        MATRIX(dxdydz, k, 0)-=xd*rf;
        MATRIX(dxdydz, j, 1)+=yd*rf;
        MATRIX(dxdydz, k, 1)-=yd*rf;
        MATRIX(dxdydz, j, 2)+=zd*rf;
        MATRIX(dxdydz, k, 2)-=zd*rf;
      }
    }

    /*Calculate the attractive "force"*/
    IGRAPH_CHECK(igraph_es_all(graph, &edgeit));
    while (!igraph_es_end(graph, &edgeit)) {
      j=igraph_es_from(graph, &edgeit);
      k=igraph_es_to(graph, &edgeit);
      xd=MATRIX(*res, j, 0)-MATRIX(*res, k, 0);
      yd=MATRIX(*res, j, 1)-MATRIX(*res, k, 1);
      zd=MATRIX(*res, j, 2)-MATRIX(*res, k, 2);
      ded=sqrt(xd*xd+yd*yd+zd*zd);  /*Get dyadic euclidean distance*/
      if (ded != 0) {
	xd/=ded;                      /*Rescale differences to length 1*/
	yd/=ded;
	zd/=ded;
      }
      af=ded*ded/frk;
      MATRIX(dxdydz, j, 0)-=xd*af;   /*Add to the position change vector*/
      MATRIX(dxdydz, k, 0)+=xd*af;
      MATRIX(dxdydz, j, 1)-=yd*af;
      MATRIX(dxdydz, k, 1)+=yd*af;
      MATRIX(dxdydz, j, 2)-=zd*af;
      MATRIX(dxdydz, k, 2)+=zd*af;
      igraph_es_next(graph, &edgeit);
    }
    igraph_es_destroy(&edgeit);

    /*Dampen motion, if needed, and move the points*/
    for(j=0;j<no_of_nodes;j++){
      ded=sqrt(MATRIX(dxdydz, j, 0)*MATRIX(dxdydz, j, 0)+
	       MATRIX(dxdydz, j, 1)*MATRIX(dxdydz, j, 1)+
	       MATRIX(dxdydz, j, 2)*MATRIX(dxdydz, j, 2));
      if(ded>t){                 /*Dampen to t*/
        ded=t/ded;
        MATRIX(dxdydz, j, 0)*=ded;
        MATRIX(dxdydz, j, 1)*=ded;
        MATRIX(dxdydz, j, 2)*=ded;
      }
      MATRIX(*res, j, 0)+=MATRIX(dxdydz, j, 0);          /*Update positions*/
      MATRIX(*res, j, 1)+=MATRIX(dxdydz, j, 1);
      MATRIX(*res, j, 2)+=MATRIX(dxdydz, j, 2);
    }
  }
  
  igraph_matrix_destroy(&dxdydz);
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

int igraph_layout_kamada_kawai(const igraph_t *graph, igraph_matrix_t *res,
			       integer_t niter, real_t sigma, 
			       real_t initemp, real_t coolexp,
			       real_t kkconst) {

  real_t temp, candx, candy;
  real_t dpot, odis, ndis, osqd, nsqd;
  long int n,i,j,k;
  igraph_matrix_t elen;

  /* Define various things */
  n=igraph_vcount(graph);

  /* Calculate elen, initial x & y */

  RNG_BEGIN();

  IGRAPH_CHECK(igraph_matrix_resize(res, n, 2));
  IGRAPH_MATRIX_INIT_FINALLY(&elen, n, n);
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
  igraph_matrix_destroy(&elen);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

/**
 * \function igraph_layout_kamada_kawai_3d
 * \brief 3D version of the force based Kamada-Kawai layout, the pair
 * of the \ref igraph_layout_kamada_kawai 2D layout generator
 * 
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
 * Added in version 0.2.</para><para>
 * 
 * Time complexity: O(|V|^2) for each
 * iteration, |V| is the number of
 * vertices in the graph. 
 */

int igraph_layout_kamada_kawai_3d(const igraph_t *graph, igraph_matrix_t *res,
				  integer_t niter, real_t sigma, 
				  real_t initemp, real_t coolexp, 
				  real_t kkconst) {
  real_t temp, candx, candy, candz;
  real_t dpot, odis, ndis, osqd, nsqd;
  long int i,j,k;
  long int no_of_nodes=igraph_vcount(graph);
  igraph_matrix_t elen;
  
  RNG_BEGIN();
  
  IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 3));
  IGRAPH_MATRIX_INIT_FINALLY(&elen,  no_of_nodes, no_of_nodes);
  IGRAPH_CHECK(igraph_shortest_paths(graph, &elen, IGRAPH_VS_ALL(graph), 
				     IGRAPH_ALL));
  
  temp=initemp;
  for(i=0;i<niter;i++){
    /*Update each vertex*/
    for(j=0;j<no_of_nodes;j++){
      /*Draw the candidate via a gaussian perturbation*/
      candx=RNG_NORMAL(MATRIX(*res, j, 0) ,sigma*temp/initemp);
      candy=RNG_NORMAL(MATRIX(*res, j, 1) ,sigma*temp/initemp);
      candz=RNG_NORMAL(MATRIX(*res, j, 2) ,sigma*temp/initemp);
      /*Calculate the potential difference for the new position*/
      dpot=0.0;
      for(k=0;k<no_of_nodes;k++) /*Potential differences for pairwise effects*/
        if(j!=k){
          odis=sqrt((MATRIX(*res, j, 0)-MATRIX(*res, k, 0))*
		    (MATRIX(*res, j, 0)-MATRIX(*res, k, 0))+
		    (MATRIX(*res, j, 1)-MATRIX(*res, k, 1))*
		    (MATRIX(*res, j, 1)-MATRIX(*res, k, 1))+
		    (MATRIX(*res, j, 2)-MATRIX(*res, k, 2))*
		    (MATRIX(*res, j, 2)-MATRIX(*res, k, 2)));
          ndis=sqrt((candx-MATRIX(*res, k, 0))*(candx-MATRIX(*res, k, 0))+
		    (candy-MATRIX(*res, k, 1))*(candy-MATRIX(*res, k, 1))+
		    (candz-MATRIX(*res, k, 2))*(candz-MATRIX(*res, k, 2)));
          osqd=(odis-MATRIX(elen, j, k))*(odis-MATRIX(elen, j, k));
          nsqd=(ndis-MATRIX(elen, j, k))*(ndis-MATRIX(elen, j, k));
          dpot+=kkconst*(osqd-nsqd)/(MATRIX(elen, j, k)*MATRIX(elen, j, k));
        }
      /*Make a keep/reject decision*/
      if(log(RNG_UNIF(0.0,1.0))<dpot/temp){
        MATRIX(*res, j, 0)=candx;
        MATRIX(*res, j, 1)=candy;
        MATRIX(*res, j, 2)=candz;
      }
    }
    /*Cool the system*/
    temp*=coolexp;
  }

  RNG_END();
  igraph_matrix_destroy(&elen);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

int igraph_layout_springs(const igraph_t *graph, igraph_matrix_t *res,
			  real_t mass, real_t equil, real_t k,
			  real_t repeqdis, real_t kfr, bool_t repulse) {
  /* TODO */
  return 0;
}

void igraph_i_norm2d(real_t *x, real_t *y) {
  real_t len=sqrt((*x)*(*x) + (*y)*(*y));
  if (len != 0) {
    *x /= len;
    *y /= len;
  }
}

/**
 * \function igraph_layout_lgl
 * \brief Force based layout algorithm for large graphs.
 * 
 * This is a layout generator similar to the Large Graph Layout
 * algorithm and program
 * (http://bioinformatics.icmb.utexas.edu/lgl/). But unlike LGL, this
 * version uses a Fruchterman-Reingold style simulated annealing
 * algorithm for placing the vertices. The speedup is achived by
 * placing the vertices on a grid and calculating the repulsion only
 * for vertices which are closer to each other than a limit. 
 * 
 * \param graph The (initialized) graph object to place.
 * \param res Pointer to an initialized matrix object to hold the
 *   result. It will be resized if needed.
 * \param maxit The maximum number of cooling iterations to perform
 *   for each layout step.
 * \param maxdelta The maximum length of the move allowed for a vertex
 *   in a single iteration. 
 * \param area This parameter gives the area of the square on which
 *   the vertices will be placed.
 * \param coolexp The cooling exponent. 
 * \param repulserad Determines the radius at which vertex-vertex 
 *   repulsion cancels out attraction of adjacenct vertices.
 * \param cellsize The size of the grid cells, one side of the
 *   square. 
 * \param proot The root vertex, this is placed first, its neighbors
 *   in the first iteration, second neighbors in the second, etc. If
 *   negative then a random vertex is chosen.
 * \return Error code.
 * 
 * Added in version 0.2.</para><para>
 * 
 * Time complexity: ideally O(dia*maxit*(|V|+|E|)), |V| is the number
 * of vertices, 
 * dia is the diameter of the graph, worst case complexity is still 
 * O(dia*maxit*(|V|^2+|E|)), this is the case when all vertices happen to be
 * in the same grid cell. 
 */

int igraph_layout_lgl(const igraph_t *graph, igraph_matrix_t *res,
		      integer_t maxit, real_t maxdelta, 
		      real_t area, real_t coolexp,
		      real_t repulserad, real_t cellsize, 
		      integer_t proot) {
  
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_t mst;
  long int root;
  long int no_of_layers, actlayer=0;
  igraph_vector_t vids;
  igraph_vector_t layers;
  igraph_vector_t parents;
  igraph_vector_t edges;
  igraph_2dgrid_t grid;  
  igraph_vector_t eids;
  igraph_vector_t forcex;
  igraph_vector_t forcey;

  real_t frk=sqrt(area/no_of_nodes);
  real_t H_n=0;

  IGRAPH_CHECK(igraph_minimum_spanning_tree_unweighted(graph, &mst));
  IGRAPH_FINALLY(igraph_destroy, &mst);
  
  IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 2));
  
  /* Determine the root vertex, random pick right now */
  if (proot < 0) {
    root=RNG_INTEGER(0, no_of_nodes-1);
  } else {
    root=proot;
  }

  /* Assign the layers */
  IGRAPH_VECTOR_INIT_FINALLY(&vids, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&layers, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&parents, 0);
  IGRAPH_CHECK(igraph_bfs(&mst, root, IGRAPH_ALL, &vids, &layers, &parents));
  no_of_layers=igraph_vector_size(&layers)-1;
  
  /* We don't need the mst any more */
  igraph_destroy(&mst);
  igraph_empty(&mst, 0, IGRAPH_UNDIRECTED); /* to make finalization work */
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_CHECK(igraph_vector_reserve(&edges, no_of_edges));
  IGRAPH_VECTOR_INIT_FINALLY(&eids, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&forcex, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&forcey, no_of_nodes);

  /* Place the vertices randomly */
  IGRAPH_CHECK(igraph_layout_random(graph, res));
  igraph_matrix_multiply(res, 1e6);
  
  /* This is the grid for calculating the vertices near to a given vertex */  
  IGRAPH_CHECK(igraph_2dgrid_init(&grid, res, 
				  -sqrt(area/M_PI),sqrt(area/M_PI), cellsize,
				  -sqrt(area/M_PI),sqrt(area/M_PI), cellsize));
  IGRAPH_FINALLY(igraph_2dgrid_destroy, &grid);

  /* Place the root vertex */
  igraph_2dgrid_add(&grid, root, 0, 0);

  for (actlayer=1; actlayer<no_of_layers; actlayer++) {
    H_n += 1.0/actlayer;
  }

  for (actlayer=1; actlayer<no_of_layers; actlayer++) {

    real_t c=1;
    long int i, j;
    real_t massx, massy;
    real_t px, py;
    real_t sx, sy;

    long int it=0;
    real_t epsilon=10e-6;
    real_t maxchange=epsilon+1;
    long int pairs;
    real_t sconst=sqrt(area/M_PI) / H_n; 
    igraph_2dgrid_iterator_t vidit;

/*     printf("Layer %li:\n", actlayer); */
    
    /*-----------------------------------------*/
    /* Step 1: place the next layer on spheres */
    /*-----------------------------------------*/

    RNG_BEGIN();

    j=VECTOR(layers)[actlayer];
    for (i=VECTOR(layers)[actlayer-1]; i<VECTOR(layers)[actlayer]; i++) {

      long int vid=VECTOR(vids)[i];
      long int par=VECTOR(parents)[vid];
      igraph_2dgrid_getcenter(&grid, &massx, &massy);
      igraph_i_norm2d(&massx, &massy);
      px=MATRIX(*res, vid, 0)-MATRIX(*res, par, 0);
      py=MATRIX(*res, vid, 1)-MATRIX(*res, par, 1);
      igraph_i_norm2d(&px, &py);
      sx=c*(massx+px)+MATRIX(*res, vid, 0);
      sy=c*(massy+py)+MATRIX(*res, vid, 1);

      /* The neighbors of 'vid' */
      while (j < VECTOR(layers)[actlayer+1] && 
	     VECTOR(parents)[(long int)VECTOR(vids)[j]]==vid) {
	real_t rx, ry;
	if (actlayer==1) {
	  real_t phi=2*M_PI/(VECTOR(layers)[2]-1)*(j-1);
	  rx=cos(phi);
	  ry=sin(phi);
	} else {
	  rx=RNG_UNIF(-1,1);
	  ry=RNG_UNIF(-1,1);
	}
	igraph_i_norm2d(&rx, &ry);
	rx = rx / actlayer * sconst;
	ry = ry / actlayer * sconst;
	igraph_2dgrid_add(&grid, VECTOR(vids)[j], sx+rx, sy+ry);
	j++; 
      }
    }

    RNG_END();
    
    /*-----------------------------------------*/
    /* Step 2: add the edges of the next layer */
    /*-----------------------------------------*/

    for (j=VECTOR(layers)[actlayer]; j<VECTOR(layers)[actlayer+1]; j++) {
      long int vid=VECTOR(vids)[j];
      long int k;
      IGRAPH_CHECK(igraph_adjacent(graph, &eids, vid, IGRAPH_ALL));
      for (k=0;k<igraph_vector_size(&eids);k++) {
	long int eid=VECTOR(eids)[k];
	integer_t from, to;
	igraph_edge(graph, eid, &from, &to);
	if (from != vid && igraph_2dgrid_in(&grid, from) ||
	    to   != vid && igraph_2dgrid_in(&grid, to)) {
	  igraph_vector_push_back(&edges, eid);
	}
      }
    }

    /*-----------------------------------------*/
    /* Step 3: let the springs spring          */
    /*-----------------------------------------*/
    
    maxchange=epsilon+1;
    while (it < maxit && maxchange > epsilon) {
      long int j;
      real_t t=maxdelta*pow((maxit-it)/(double)maxit, coolexp);
      long int vid, nei;

      /* init */
      igraph_vector_null(&forcex);
      igraph_vector_null(&forcey);
      maxchange=0;
      
      /* attractive "forces" along the edges */
      for (j=0; j<igraph_vector_size(&edges); j++) {
	integer_t from, to;
	real_t xd, yd, dist, force;
	igraph_edge(graph, VECTOR(edges)[j], &from, &to);
	xd=MATRIX(*res, (long int)from, 0)-MATRIX(*res, (long int)to, 0);
	yd=MATRIX(*res, (long int)from, 1)-MATRIX(*res, (long int)to, 1);
	dist=sqrt(xd*xd+yd*yd);
	if (dist != 0) { xd/=dist; yd/=dist; }
	force=dist*dist/frk;
	VECTOR(forcex)[(long int)from] -= xd*force;
	VECTOR(forcex)[(long int)to]   += xd*force;
	VECTOR(forcey)[(long int)from] -= yd*force;
	VECTOR(forcey)[(long int)to]   += yd*force;
      }
      
      /* repulsive "forces" of the vertices nearby */
      pairs=0;
      igraph_2dgrid_reset(&grid, &vidit);
      while ( (vid=igraph_2dgrid_next(&grid, &vidit)-1) != -1) {
	while ( (nei=igraph_2dgrid_next_nei(&grid, &vidit)-1) != -1) {
	  real_t xd=MATRIX(*res, (long int)vid, 0)-
	    MATRIX(*res, (long int)nei, 0);
	  real_t yd=MATRIX(*res, (long int)vid, 1)-
	    MATRIX(*res, (long int)nei, 1);
	  real_t dist=sqrt(xd*xd+yd*yd);
	  real_t force;
	  if (dist < cellsize) {
	    pairs++;
	    if (dist==0) { dist=epsilon; };
	    xd/=dist; yd/=dist;
	    force=frk*frk*(1.0/dist-dist*dist/repulserad);
	    VECTOR(forcex)[(long int)vid] += xd*force;
	    VECTOR(forcex)[(long int)nei] -= xd*force;
	    VECTOR(forcey)[(long int)vid] += yd*force;
	    VECTOR(forcey)[(long int)nei] -= yd*force;	    
	  }
	}
      }
      
/*       printf("verties: %li iterations: %li\n",  */
/* 	     (long int) VECTOR(layers)[actlayer+1], pairs); */
      
      /* apply the changes */
      for (j=0; j<VECTOR(layers)[actlayer+1]; j++) {
	long int vid=VECTOR(vids)[j];
	real_t fx=VECTOR(forcex)[vid];
	real_t fy=VECTOR(forcey)[vid];
	real_t ded=sqrt(fx*fx+fy*fy);
	if (ded > t) {
	  ded=t/ded;
	  fx*=ded; fy *=ded;
	}
	igraph_2dgrid_move(&grid, vid, fx, fy);
	if (fx > maxchange) { maxchange=fx; }
	if (fy > maxchange) { maxchange=fy; }
      }
      it++;
/*       printf("%li iterations, maxchange: %f\n", it, (double)maxchange); */
    }
  }

  igraph_destroy(&mst);
  igraph_vector_destroy(&vids);
  igraph_vector_destroy(&layers);
  igraph_vector_destroy(&parents);
  igraph_vector_destroy(&edges);
  igraph_2dgrid_destroy(&grid);
  igraph_vector_destroy(&eids);
  igraph_vector_destroy(&forcex);
  igraph_vector_destroy(&forcey);
  IGRAPH_FINALLY_CLEAN(9);
  return 0;

}

/**
 * \function igraph_layout_grid_fruchterman_reingold
 * \brief Force based layout generator for large graphs.
 * 
 * This algorithm is the same as the Fruchterman-Reingold layout
 * generator, but it partitions the 2d space to a grid and and vertex
 * repulsion is calculated only for vertices nearby.
 *
 * \param graph The graph object. 
 * \param res The result, the coordinates in a matrix. The parameter
 *   should point to an initialized matrix object and will be resized.
 * \param niter Number of iterations to perform.
 * \param maxdelta Maximum distance to move a vertex in an iteration.
 * \param area The area of the square on which the vertices will be
 *   placed.
 * \param coolexp The cooling exponent.
 * \param repulserad Determines the radius at which vertex-vertex 
 *   repulsion cancels out attraction of adjacenct vertices.
 * \param cellsize The size of the grid cells.
 * \param use_seed Logical, if true, the coordinates passed in \p res
 *   (should have the appropriate size) will be used for the first
 *   iteration.
 * \return Error code.
 *
 * Added in version 0.2.</para><para>
 * 
 * Time complexity: ideally (constant number of vertices in each cell) 
 * O(niter*(|V|+|E|)), in the worst case O(niter*(|V|^2+|E|)).
 */

int igraph_layout_grid_fruchterman_reingold(const igraph_t *graph, 
					    igraph_matrix_t *res,
					    integer_t niter, real_t maxdelta, 
					    real_t area, real_t coolexp,
					    real_t repulserad, 
					    real_t cellsize,
					    bool_t use_seed) {

  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_2dgrid_t grid;  
  igraph_vector_t forcex;
  igraph_vector_t forcey;
  long int i, it=0;
  igraph_2dgrid_iterator_t vidit;  

  real_t frk=sqrt(area/no_of_nodes);

  IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 2));
  IGRAPH_VECTOR_INIT_FINALLY(&forcex, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&forcey, no_of_nodes);
  
  /* initial layout */
  if (!use_seed) {
    IGRAPH_CHECK(igraph_layout_random(graph, res));
    igraph_matrix_multiply(res, sqrt(area/M_PI));
  }
  
  /* make grid */
  IGRAPH_CHECK(igraph_2dgrid_init(&grid, res, 
				  -sqrt(area/M_PI),sqrt(area/M_PI), cellsize,
				  -sqrt(area/M_PI),sqrt(area/M_PI), cellsize));
  IGRAPH_FINALLY(igraph_2dgrid_destroy, &grid);
  
  /* place vertices on grid */
  for (i=0; i<no_of_nodes; i++) {
    igraph_2dgrid_add2(&grid, i);
  }

  while (it<niter) {
    long int j;
    real_t t=maxdelta*pow((niter-it)/(double)niter, coolexp);
    long int vid, nei;

    igraph_vector_null(&forcex);
    igraph_vector_null(&forcey);
    
    /* attraction */
    for (j=0; j<no_of_edges; j++) {
      integer_t from, to;
      real_t xd, yd, dist, force;
      igraph_edge(graph, j, &from, &to);
      xd=MATRIX(*res, (long int)from, 0)-MATRIX(*res, (long int)to, 0);
      yd=MATRIX(*res, (long int)from, 1)-MATRIX(*res, (long int)to, 1);
      dist=sqrt(xd*xd+yd*yd);
      if (dist != 0) { xd/=dist; yd/=dist; }
      force=dist*dist/frk;
      VECTOR(forcex)[(long int)from] -= xd*force;
      VECTOR(forcex)[(long int)to]   += xd*force;
      VECTOR(forcey)[(long int)from] -= yd*force;
      VECTOR(forcey)[(long int)to]   += yd*force;
    }

    /* repulsion */
    igraph_2dgrid_reset(&grid, &vidit);
    while ( (vid=igraph_2dgrid_next(&grid, &vidit)-1) != -1) {
      while ( (nei=igraph_2dgrid_next_nei(&grid, &vidit)-1) != -1) {
	real_t xd=MATRIX(*res, (long int)vid, 0)-
	  MATRIX(*res, (long int)nei, 0);
	real_t yd=MATRIX(*res, (long int)vid, 1)-
	  MATRIX(*res, (long int)nei, 1);
	real_t dist=sqrt(xd*xd+yd*yd);
	real_t force;
	if (dist < cellsize) {
	  if (dist==0) { dist=1e-6; };
	  xd/=dist; yd/=dist;
	  force=frk*frk*(1.0/dist-dist*dist/repulserad);
	  VECTOR(forcex)[(long int)vid] += xd*force;
	  VECTOR(forcex)[(long int)nei] -= xd*force;
	  VECTOR(forcey)[(long int)vid] += yd*force;
	  VECTOR(forcey)[(long int)nei] -= yd*force;	    
	}
      }
    }

    /* update */
    for (j=0; j<no_of_nodes; j++) {
      long int vid=j;
      real_t fx=VECTOR(forcex)[vid];
      real_t fy=VECTOR(forcey)[vid];
      real_t ded=sqrt(fx*fx+fy*fy);
      if (ded > t) {
	ded=t/ded;
	fx*=ded; fy *=ded;
      }
      igraph_2dgrid_move(&grid, vid, fx, fy);
    }
    it++;
    
  } /* it<niter */

  IGRAPH_FINALLY_CLEAN(3);
  return 0;
}
