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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include "igraph.h"
#include "random.h"
#include "memory.h"
#include <math.h>

int igraph_i_layout_sphere_2d(igraph_matrix_t *coords, igraph_real_t *x, igraph_real_t *y,
			      igraph_real_t *r);
int igraph_i_layout_sphere_3d(igraph_matrix_t *coords, igraph_real_t *x, igraph_real_t *y,
			      igraph_real_t *z, igraph_real_t *r);

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
    igraph_real_t phi=2*M_PI/no_of_nodes*i;
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
  igraph_real_t h;
  
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
    igraph_real_t x=cos(MATRIX(*res, i, 1))*sin(MATRIX(*res, i, 0));
    igraph_real_t y=sin(MATRIX(*res, i, 1))*sin(MATRIX(*res, i, 0));
    igraph_real_t z=cos(MATRIX(*res, i, 0));
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
				       igraph_integer_t niter, igraph_real_t maxdelta,
				       igraph_real_t area, igraph_real_t coolexp, 
				       igraph_real_t repulserad, igraph_bool_t use_seed) {
  igraph_real_t frk,t,ded,xd,yd;
  igraph_real_t rf,af;
  long int i,j,k;

  long int no_of_nodes=igraph_vcount(graph);
  igraph_matrix_t dxdy=IGRAPH_MATRIX_NULL;
  igraph_eit_t edgeit;
  igraph_integer_t from, to;
  
  IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 2));
  if (!use_seed) {
    IGRAPH_CHECK(igraph_layout_random(graph, res));
  }
  IGRAPH_MATRIX_INIT_FINALLY(&dxdy, no_of_nodes, 2);

  IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(0), &edgeit));
  IGRAPH_FINALLY(igraph_eit_destroy, &edgeit);
  
  frk=sqrt(area/no_of_nodes);

  for(i=niter;i>0;i--) {
    /* Report progress in approx. every 100th step */
    if (i%10 == 0)
      igraph_progress("Fruchterman-Reingold layout: ",
		      100.0-100.0*i/niter, NULL);
    
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
    IGRAPH_EIT_RESET(edgeit);
    while (!IGRAPH_EIT_END(edgeit)) {
      igraph_edge(graph, IGRAPH_EIT_GET(edgeit), &from, &to);
      j=from; 
      k=to;
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
      IGRAPH_EIT_NEXT(edgeit);
    }
    
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

  igraph_progress("Fruchterman-Reingold layout: ", 100.0, NULL);
  
  igraph_eit_destroy(&edgeit);
  igraph_matrix_destroy(&dxdy);
  IGRAPH_FINALLY_CLEAN(2);
  
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
					  igraph_integer_t niter, igraph_real_t maxdelta,
					  igraph_real_t volume, igraph_real_t coolexp,
					  igraph_real_t repulserad,
					  igraph_bool_t use_seed) {
  
  igraph_real_t frk, t, ded, xd, yd, zd;
  igraph_matrix_t dxdydz;
  igraph_real_t rf, af;
  long int i, j, k;
  
  long int no_of_nodes=igraph_vcount(graph);
  igraph_eit_t edgeit;
  igraph_integer_t from, to;
  
  IGRAPH_CHECK(igraph_matrix_init(&dxdydz, no_of_nodes, 3));
  IGRAPH_FINALLY(igraph_matrix_destroy, &dxdydz);
  
  IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 3));
  if (!use_seed) {
    IGRAPH_CHECK(igraph_layout_random_3d(graph, res));
  }

  IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(0), &edgeit));
  IGRAPH_FINALLY(igraph_eit_destroy, &edgeit);

  frk=pow(volume/(double)no_of_nodes,1.0/3.0); /*Define the F-R constant*/

  /*Run the annealing loop*/
  for(i=niter;i>=0;i--){
    if (i%10 == 0)
      igraph_progress("3D Fruchterman-Reingold layout: ",
		      100.0-100.0*i/niter, NULL);

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
    IGRAPH_EIT_RESET(edgeit);
    while (!IGRAPH_EIT_END(edgeit)) {
      igraph_edge(graph, IGRAPH_EIT_GET(edgeit), &from, &to);
      j=from;
      k=to;
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
      IGRAPH_EIT_NEXT(edgeit);
    }

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

  igraph_progress("3D Fruchterman-Reingold layout: ", 100.0, NULL);
  
  igraph_matrix_destroy(&dxdydz);
  igraph_eit_destroy(&edgeit);
  IGRAPH_FINALLY_CLEAN(2);
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
			       igraph_integer_t niter, igraph_real_t sigma, 
			       igraph_real_t initemp, igraph_real_t coolexp,
			       igraph_real_t kkconst) {

  igraph_real_t temp, candx, candy;
  igraph_real_t dpot, odis, ndis, osqd, nsqd;
  long int n,i,j,k;
  igraph_matrix_t elen;

  /* Define various things */
  n=igraph_vcount(graph);

  /* Calculate elen, initial x & y */

  RNG_BEGIN();

  IGRAPH_CHECK(igraph_matrix_resize(res, n, 2));
  IGRAPH_MATRIX_INIT_FINALLY(&elen, n, n);
  IGRAPH_CHECK(igraph_shortest_paths(graph, &elen, igraph_vss_all(), 
				     IGRAPH_ALL));
  for (i=0; i<n; i++) {
    MATRIX(elen, i, i) = sqrt(n);
    MATRIX(*res, i, 0) = RNG_NORMAL(0, n/4.0);
    MATRIX(*res, i, 1) = RNG_NORMAL(0, n/4.0);
  }
  
  /*Perform the annealing loop*/
  temp=initemp;
  for(i=0;i<niter;i++){
    /* Report progress in approx. every 100th step */
    if (i%10 == 0)
      igraph_progress("Kamada-Kawai layout: ",
		      100.0*i/niter, NULL);
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

  igraph_progress("Kamada-Kawai layout: ", 100.0, NULL);

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
				  igraph_integer_t niter, igraph_real_t sigma, 
				  igraph_real_t initemp, igraph_real_t coolexp, 
				  igraph_real_t kkconst) {
  igraph_real_t temp, candx, candy, candz;
  igraph_real_t dpot, odis, ndis, osqd, nsqd;
  long int i,j,k;
  long int no_of_nodes=igraph_vcount(graph);
  igraph_matrix_t elen;
  
  RNG_BEGIN();
  
  IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 3));
  IGRAPH_MATRIX_INIT_FINALLY(&elen,  no_of_nodes, no_of_nodes);
  IGRAPH_CHECK(igraph_shortest_paths(graph, &elen, igraph_vss_all(), 
				     IGRAPH_ALL));
  
  temp=initemp;
  for(i=0;i<niter;i++){
    if (i%10 == 0)
      igraph_progress("3D Kamada-Kawai layout: ",
		      100.0*i/niter, NULL);

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

  igraph_progress("3D Kamada-Kawai layout: ", 100.0, NULL);

  RNG_END();
  igraph_matrix_destroy(&elen);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

int igraph_layout_springs(const igraph_t *graph, igraph_matrix_t *res,
			  igraph_real_t mass, igraph_real_t equil, igraph_real_t k,
			  igraph_real_t repeqdis, igraph_real_t kfr, igraph_bool_t repulse) {
  /* TODO */
  return 0;
}

void igraph_i_norm2d(igraph_real_t *x, igraph_real_t *y) {
  igraph_real_t len=sqrt((*x)*(*x) + (*y)*(*y));
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
		      igraph_integer_t maxit, igraph_real_t maxdelta, 
		      igraph_real_t area, igraph_real_t coolexp,
		      igraph_real_t repulserad, igraph_real_t cellsize, 
		      igraph_integer_t proot) {
  
  
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

  igraph_real_t frk=sqrt(area/no_of_nodes);
  igraph_real_t H_n=0;

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

    igraph_real_t c=1;
    long int i, j;
    igraph_real_t massx, massy;
    igraph_real_t px, py;
    igraph_real_t sx, sy;

    long int it=0;
    igraph_real_t epsilon=10e-6;
    igraph_real_t maxchange=epsilon+1;
    long int pairs;
    igraph_real_t sconst=sqrt(area/M_PI) / H_n; 
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
	igraph_real_t rx, ry;
	if (actlayer==1) {
	  igraph_real_t phi=2*M_PI/(VECTOR(layers)[2]-1)*(j-1);
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
	igraph_integer_t from, to;
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
      igraph_real_t t=maxdelta*pow((maxit-it)/(double)maxit, coolexp);
      long int vid, nei;

      /* init */
      igraph_vector_null(&forcex);
      igraph_vector_null(&forcey);
      maxchange=0;
      
      /* attractive "forces" along the edges */
      for (j=0; j<igraph_vector_size(&edges); j++) {
	igraph_integer_t from, to;
	igraph_real_t xd, yd, dist, force;
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
	  igraph_real_t xd=MATRIX(*res, (long int)vid, 0)-
	    MATRIX(*res, (long int)nei, 0);
	  igraph_real_t yd=MATRIX(*res, (long int)vid, 1)-
	    MATRIX(*res, (long int)nei, 1);
	  igraph_real_t dist=sqrt(xd*xd+yd*yd);
	  igraph_real_t force;
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
	igraph_real_t fx=VECTOR(forcex)[vid];
	igraph_real_t fy=VECTOR(forcey)[vid];
	igraph_real_t ded=sqrt(fx*fx+fy*fy);
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
					    igraph_integer_t niter, igraph_real_t maxdelta, 
					    igraph_real_t area, igraph_real_t coolexp,
					    igraph_real_t repulserad, 
					    igraph_real_t cellsize,
					    igraph_bool_t use_seed) {

  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_2dgrid_t grid;  
  igraph_vector_t forcex;
  igraph_vector_t forcey;
  long int i, it=0;
  igraph_2dgrid_iterator_t vidit;  

  igraph_real_t frk=sqrt(area/no_of_nodes);

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
    igraph_real_t t=maxdelta*pow((niter-it)/(double)niter, coolexp);
    long int vid, nei;

    /* Report progress */
    if (it%10 == 0) {
      igraph_progress("Grid based Fruchterman-Reingold layout: ", 
		      (100.0*it)/niter, NULL);
    }

    igraph_vector_null(&forcex);
    igraph_vector_null(&forcey);
    
    /* attraction */
    for (j=0; j<no_of_edges; j++) {
      igraph_integer_t from, to;
      igraph_real_t xd, yd, dist, force;
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
	igraph_real_t xd=MATRIX(*res, (long int)vid, 0)-
	  MATRIX(*res, (long int)nei, 0);
	igraph_real_t yd=MATRIX(*res, (long int)vid, 1)-
	  MATRIX(*res, (long int)nei, 1);
	igraph_real_t dist=sqrt(xd*xd+yd*yd);
	igraph_real_t force;
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
      igraph_real_t fx=VECTOR(forcex)[vid];
      igraph_real_t fy=VECTOR(forcey)[vid];
      igraph_real_t ded=sqrt(fx*fx+fy*fy);
      if (ded > t) {
	ded=t/ded;
	fx*=ded; fy *=ded;
      }
      igraph_2dgrid_move(&grid, vid, fx, fy);
    }
    it++;
    
  } /* it<niter */
  
  igraph_progress("Grid based Fruchterman-Reingold layout: ", 
		  100.0, NULL);

  igraph_vector_destroy(&forcex);
  igraph_vector_destroy(&forcey);
  igraph_2dgrid_destroy(&grid);
  IGRAPH_FINALLY_CLEAN(3);
  return 0;
}

/* Internal structure for Reingold-Tilford layout */
struct igraph_i_reingold_tilford_vertex {
  int parent;        /* Parent node index */
  int level;         /* Level of the node */
  igraph_real_t offset;     /* X offset from parent node */
  int left_contour;  /* Next left node of the contour
		      of the subtree rooted at this node */
  int right_contour; /* Next right node of the contour
		      of the subtree rooted at this node */
};

int igraph_i_layout_reingold_tilford_postorder(struct igraph_i_reingold_tilford_vertex *vdata,
                                               long int node, long int vcount);
int igraph_i_layout_reingold_tilford_calc_coords(struct igraph_i_reingold_tilford_vertex *vdata,
                                                 igraph_matrix_t *res, long int node,
												 long int vcount, igraph_real_t xpos);

/**
 * \function igraph_layout_reingold_tilford
 * \brief Reingold-Tilford layout for tree graphs
 * 
 * Arranges the nodes in a tree where the given node is used as the root.
 * The tree is directed downwards and the parents are centered above its
 * children. For the exact algorithm, see:
 * 
 * Reingold, E and Tilford, J: Tidier drawing of trees.
 * IEEE Trans. Softw. Eng., SE-7(2):223--228, 1981
 *
 * If the given graph is not a tree, a breadth-first search is executed
 * first to obtain a possible spanning tree.
 * 
 * \param graph The graph object. 
 * \param res The result, the coordinates in a matrix. The parameter
 *   should point to an initialized matrix object and will be resized.
 * \param root The index of the root vertex.
 * \return Error code.
 *
 * Added in version 0.2.
 * 
 * TODO: decompose and merge for not fully connected graphs
 * TODO: possible speedup could be achieved if we use a table for storing
 * the children of each node in the tree. (Now the implementation uses a
 * single array containing the parent of each node and a node's children
 * are determined by looking for other nodes that have this node as parent)
 */
int igraph_layout_reingold_tilford(const igraph_t *graph, 
				   igraph_matrix_t *res, long int root) {
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  long int i, n, j, it=0;
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  igraph_i_adjlist_t allneis;
  igraph_vector_t *neis;
  struct igraph_i_reingold_tilford_vertex *vdata;
  
  if (root<0 || root>=no_of_nodes) {
    IGRAPH_ERROR("invalid vertex id", IGRAPH_EINVVID);
  }
  
  IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 2));
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
  
  IGRAPH_CHECK(igraph_i_adjlist_init(graph, &allneis, IGRAPH_ALL));
  IGRAPH_FINALLY(igraph_i_adjlist_destroy, &allneis);
  
  vdata=Calloc(no_of_nodes, struct igraph_i_reingold_tilford_vertex);
  if (vdata==0) {
    IGRAPH_ERROR("igraph_layout_reingold_tilford failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, vdata);

  for (i=0; i<no_of_nodes; i++) {
    vdata[i].parent=-1;
    vdata[i].level=-1;
    vdata[i].offset=0.0;
    vdata[i].left_contour=-1;
    vdata[i].right_contour=-1;
  }
  vdata[root].parent=root;
  vdata[root].level=0;
  MATRIX(*res, root, 1) = 0;
  
  /* Step 1: assign Y coordinates based on BFS and setup parents vector */
  IGRAPH_CHECK(igraph_dqueue_push(&q, root));
  IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
  while (!igraph_dqueue_empty(&q)) {
    long int actnode=igraph_dqueue_pop(&q);
    long int actdist=igraph_dqueue_pop(&q);
    neis=igraph_i_adjlist_get(&allneis, actnode);
    n=igraph_vector_size(neis);
    for (j=0; j<n; j++) {
      long int neighbor=VECTOR(*neis)[j];
      if (vdata[neighbor].parent >= 0) { continue; }
      MATRIX(*res, neighbor, 1)=actdist+1;
      IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
      IGRAPH_CHECK(igraph_dqueue_push(&q, actdist+1));
      vdata[neighbor].parent = actnode;
      vdata[neighbor].level = actdist+1;
    }
  }
  
  /* Step 2: postorder tree traversal, determines the appropriate X
   * offsets for every node */
  igraph_i_layout_reingold_tilford_postorder(vdata, root, no_of_nodes);
  
  /* Step 3: calculate real coordinates based on X offsets */
  igraph_i_layout_reingold_tilford_calc_coords(vdata, res, root, no_of_nodes, vdata[root].offset);
  
  igraph_dqueue_destroy(&q);
  igraph_i_adjlist_destroy(&allneis);
  igraph_free(vdata);
  IGRAPH_FINALLY_CLEAN(3);
  
  igraph_progress("Reingold-Tilford tree layout", 100.0, NULL);
  
  return 0;
}

int igraph_i_layout_reingold_tilford_calc_coords(struct igraph_i_reingold_tilford_vertex *vdata,
                                                 igraph_matrix_t *res, long int node,
												 long int vcount, igraph_real_t xpos) {
  long int i, n;
  igraph_real_t avg=0.0;
  MATRIX(*res, node, 0) = xpos;
  for (i=0, n=0; i<vcount; i++) {
    if (i == node) continue;
    if (vdata[i].parent == node) {
      igraph_i_layout_reingold_tilford_calc_coords(vdata, res, i, vcount,
						   xpos+vdata[i].offset);
    }
  }
}

int igraph_i_layout_reingold_tilford_postorder(struct igraph_i_reingold_tilford_vertex *vdata,
                                               long int node, long int vcount) {
  long int i, j, childcount, leftroot, leftrootidx;
  igraph_real_t avg;
  
  /* printf("Starting visiting node %d\n", node); */
  
  /* Check whether this node is a leaf node */
  childcount=0;
  for (i=0; i<vcount; i++) {
    if (i == node) continue;
    if (vdata[i].parent == node) {
      /* Node i is a child, so visit it recursively */
      childcount++;
      igraph_i_layout_reingold_tilford_postorder(vdata, i, vcount);
    }
  }
  
  if (childcount == 0) return 0;
  
  /* Here we can assume that all of the subtrees have been placed and their
   * left and right contours are calculated. Let's place them next to each
   * other as close as we can.
   * We will take each subtree in an arbitrary order. The root of the
   * first one will be placed at offset 0, the next ones will be placed
   * as close to each other as possible. leftroot stores the root of the
   * rightmost subtree of the already placed subtrees - its right contour
   * will be checked against the left contour of the next subtree */
  leftroot=leftrootidx=-1;
  avg=0.0;
  /* printf("Visited node %d and arranged its subtrees\n", node); */
  for (i=0, j=0; i<vcount; i++) {
    if (i == node) continue;
    if (vdata[i].parent == node) {
      /* printf("  Placing child %d on level %d\n", i, vdata[i].level); */
      if (leftroot >= 0) {
	/* Now we will follow the right contour of leftroot and the
	 * left contour of the subtree rooted at i */
	long lnode, rnode;
	igraph_real_t loffset, roffset, minsep, rootsep;
	lnode = leftroot; rnode = i;
	minsep = 1;
	rootsep = vdata[leftroot].offset + minsep;
	loffset = 0; roffset = minsep;
	/* printf("    Contour: [%d, %d], offsets: [%lf, %lf], rootsep: %lf\n",
	 lnode, rnode, loffset, roffset, rootsep);*/
	while ((lnode >= 0) && (rnode >= 0)) {
	  /* Step to the next level on the left contour */
	  if (vdata[lnode].right_contour >= 0) {
	    lnode = vdata[lnode].right_contour;
	    loffset += vdata[lnode].offset;
	  } else {
	    lnode = vdata[lnode].left_contour;
	    if (lnode >= 0) loffset += vdata[lnode].offset;
	  }
	  /* Step to the next level on the right contour */
	  if (vdata[rnode].left_contour >= 0) {
	    rnode = vdata[rnode].left_contour;
	    roffset += vdata[rnode].offset;
	  } else if (vdata[rnode].right_contour >= 0) {
	    rnode = vdata[rnode].right_contour;
	    roffset += vdata[rnode].offset;
	  } else {
	    /* Right subtree ended. If the left subtree is deeper, the
	     * left contour will continue on the left subtree. We can
	     * use only lnode here, since it has already been advanced
	     * to the next level in the previous if statement */
	    vdata[rnode].left_contour = lnode;
	    vdata[rnode].right_contour = lnode;
	    rnode = -1;
	    /* printf("      Right subtree ended, continuing its contours to %d\n", vdata[rnode].left_contour); */
	  }
	  /* printf("    Contour: [%d, %d], offsets: [%lf, %lf], rootsep: %lf\n", 
	   lnode, rnode, loffset, roffset, rootsep);*/
	  
	  /* Push subtrees away if necessary */
	  if ((lnode >= 0) && (rnode >= 0) && (roffset - loffset < minsep)) {
	    /* printf("    Pushing right subtree away by %lf\n", minsep-roffset+loffset); */
	    rootsep += minsep-roffset+loffset;
	    roffset = loffset+minsep;
	  }
	}
	/* printf("  Offset of subtree with root node %d will be %lf\n", i, rootsep); */
	vdata[i].offset = rootsep;
	vdata[node].right_contour=i;
	avg = (avg*j)/(j+1) + rootsep/(j+1);
	leftrootidx=j;
	leftroot=i;
      } else {
	leftrootidx=j;
	leftroot=i;
	vdata[node].left_contour=i;
	avg = vdata[i].offset;
      }
      j++;
    }
  }
  /* printf("Shifting node to be centered above children. Shift amount: %lf\n", avg); */
  vdata[node].offset += avg;
  for (i=0, j=0; i<vcount; i++) {
    if (i == node) continue;
    if (vdata[i].parent == node) vdata[i].offset -= avg;
  }
  
  return 0;
}

int igraph_i_layout_merge_dla(igraph_i_layout_mergegrid_t *grid, 
			      long int actg, igraph_real_t *x, igraph_real_t *y, igraph_real_t r,
			      igraph_real_t cx, igraph_real_t cy, igraph_real_t startr, 
			      igraph_real_t killr);

int igraph_layout_merge_dla(igraph_vector_ptr_t *thegraphs,
			    igraph_vector_ptr_t *coords, 
			    igraph_matrix_t *res) {
  long int graphs=igraph_vector_ptr_size(coords);
  igraph_vector_t sizes;
  igraph_vector_t x, y, r;
  igraph_vector_t nx, ny, nr;
  long int allnodes=0;
  long int i, j;
  long int actg;
  igraph_i_layout_mergegrid_t grid;
  long int jpos=0;
  igraph_real_t minx, maxx, miny, maxy;
  igraph_real_t area=0;
  igraph_real_t maxr=0;
  long int respos;
  
  IGRAPH_VECTOR_INIT_FINALLY(&sizes, graphs);
  IGRAPH_VECTOR_INIT_FINALLY(&x, graphs);
  IGRAPH_VECTOR_INIT_FINALLY(&y, graphs);
  IGRAPH_VECTOR_INIT_FINALLY(&r, graphs);
  IGRAPH_VECTOR_INIT_FINALLY(&nx, graphs);
  IGRAPH_VECTOR_INIT_FINALLY(&ny, graphs);
  IGRAPH_VECTOR_INIT_FINALLY(&nr, graphs);
  
  for (i=0; i<igraph_vector_ptr_size(coords); i++) {
    igraph_matrix_t *mat=VECTOR(*coords)[i];
    long int size=igraph_matrix_nrow(mat);
    allnodes += size;
    VECTOR(sizes)[i]=size;
    VECTOR(r)[i]=pow(size, .75);
    area+=VECTOR(r)[i] * VECTOR(r)[i];
    if (VECTOR(r)[i] > maxr) {
      maxr=VECTOR(r)[i];
    }

    igraph_i_layout_sphere_2d(mat,
			      igraph_vector_e_ptr(&nx, i),
			      igraph_vector_e_ptr(&ny, i),
			      igraph_vector_e_ptr(&nr, i));
    
  }
  igraph_vector_order2(&sizes);	/* largest first */

  /* 0. create grid */
  minx=miny=-sqrt(5*area);
  maxx=maxy=sqrt(5*area);
  igraph_i_layout_mergegrid_init(&grid, minx, maxx, 200,
				 miny, maxy, 200);
  IGRAPH_FINALLY(igraph_i_layout_mergegrid_destroy, &grid);

/*   fprintf(stderr, "Ok, starting DLA\n"); */
  
  /* 1. place the largest  */
  actg=VECTOR(sizes)[jpos++];
  igraph_i_layout_merge_place_sphere(&grid, 0, 0, VECTOR(r)[actg], actg);
  
  igraph_progress("Merging layouts via DLA", 0.0, NULL);
  while (jpos<graphs) {
/*     fprintf(stderr, "comp: %li", jpos); */
    igraph_progress("Merging layouts via DLA", (100.0*jpos)/graphs, NULL);
    
    actg=VECTOR(sizes)[jpos++];
    /* 2. random walk, TODO: tune parameters */
    igraph_i_layout_merge_dla(&grid, actg, 
			      igraph_vector_e_ptr(&x, actg),
			      igraph_vector_e_ptr(&y, actg), 
			      VECTOR(r)[actg], 0, 0,
			      maxx-maxr, maxx-maxr+5);
    
    /* 3. place sphere */
    igraph_i_layout_merge_place_sphere(&grid, VECTOR(x)[actg], VECTOR(y)[actg],
				       VECTOR(r)[actg], actg);
  }
  igraph_progress("Merging layouts via DLA", 100.0, NULL);

  /* Create the result */
  IGRAPH_CHECK(igraph_matrix_resize(res, allnodes, 2));
  respos=0;
  for (i=0; i<graphs; i++) {
    long int size=igraph_matrix_nrow(VECTOR(*coords)[i]);
    igraph_real_t xx=VECTOR(x)[i];
    igraph_real_t yy=VECTOR(y)[i];
    igraph_real_t rr=VECTOR(r)[i]/VECTOR(nr)[i];
    igraph_matrix_t *mat=VECTOR(*coords)[i];
    if (VECTOR(nr)[i]==0) { rr=1; }
    for (j=0; j<size; j++) {
      MATRIX(*res, respos, 0)=rr*(MATRIX(*mat, j, 0)-VECTOR(nx)[i]);
      MATRIX(*res, respos, 1)=rr*(MATRIX(*mat, j, 1)-VECTOR(ny)[i]);
      MATRIX(*res, respos, 0)+=xx;
      MATRIX(*res, respos, 1)+=yy;
      ++respos;
    }
  }
 
  igraph_i_layout_mergegrid_destroy(&grid);
  igraph_vector_destroy(&sizes);
  igraph_vector_destroy(&x);
  igraph_vector_destroy(&y);
  igraph_vector_destroy(&r);
  igraph_vector_destroy(&nx);
  igraph_vector_destroy(&ny);
  igraph_vector_destroy(&nr);
  IGRAPH_FINALLY_CLEAN(8);
  return 0;
}

int igraph_i_layout_sphere_2d(igraph_matrix_t *coords, igraph_real_t *x, igraph_real_t *y,
			      igraph_real_t *r) {
  long int nodes=igraph_matrix_nrow(coords);
  long int i;
  igraph_real_t xmin, xmax, ymin, ymax;
  
  xmin=xmax=MATRIX(*coords,0,0);
  ymin=ymax=MATRIX(*coords,0,1);
  for (i=1; i<nodes; i++) {

    if (MATRIX(*coords,i,0) < xmin) {
      xmin=MATRIX(*coords,i,0);
    } else if (MATRIX(*coords,i,0)>xmax) {
      xmax=MATRIX(*coords,i,0);
    }

    if (MATRIX(*coords,i,1) < ymin) {
      ymin=MATRIX(*coords,i,1);
    } else if (MATRIX(*coords,i,1)>ymax) {
      ymax=MATRIX(*coords,i,1);
    }
    
  }

  *x=(xmin+xmax)/2;
  *y=(ymin+ymax)/2;
  *r=sqrt( (xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin) ) / 2;

  return 0;
}

int igraph_i_layout_sphere_3d(igraph_matrix_t *coords, igraph_real_t *x, igraph_real_t *y,
			      igraph_real_t *z, igraph_real_t *r) {
  long int nodes=igraph_matrix_nrow(coords);
  long int i;
  igraph_real_t xmin, xmax, ymin, ymax, zmin, zmax;
  
  xmin=xmax=MATRIX(*coords,0,0);
  ymin=ymax=MATRIX(*coords,0,1);
  zmin=zmax=MATRIX(*coords,0,2);
  for (i=1; i<nodes; i++) {
    
    if (MATRIX(*coords,i,0) < xmin) {
      xmin=MATRIX(*coords,i,0);
    } else if (MATRIX(*coords,i,0)>xmax) {
      xmax=MATRIX(*coords,i,0);
    }

    if (MATRIX(*coords,i,1) < ymin) {
      ymin=MATRIX(*coords,i,1);
    } else if (MATRIX(*coords,i,1)>ymax) {
      ymax=MATRIX(*coords,i,1);
    }
    
    if (MATRIX(*coords,i,2) < zmin) {
      zmin=MATRIX(*coords,i,2);
    } else if (MATRIX(*coords,i,2)>zmax) {
      zmax=MATRIX(*coords,i,2);
    }

  }
  
  *x=(xmin+xmax)/2;
  *y=(ymin+ymax)/2;
  *z=(zmin+zmax)/2;
  *r=sqrt( (xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin)+
	   (zmax-zmin)*(zmax-zmin) ) / 2;
  
  return 0;
}

#define DIST(x,y) (sqrt(pow((x)-cx,2)+pow((y)-cy,2)))

int igraph_i_layout_merge_dla(igraph_i_layout_mergegrid_t *grid, 
			      long int actg, igraph_real_t *x, igraph_real_t *y, igraph_real_t r,
			      igraph_real_t cx, igraph_real_t cy, igraph_real_t startr, 
			      igraph_real_t killr) {
  long int sp=-1;
  igraph_real_t angle, len;
  long int steps=0;

  while (sp < 0) {
    /* start particle */
    do {
      steps++;
      angle=RNG_UNIF(0,2*M_PI);
      len=RNG_UNIF(.5*startr, startr);
      *x=cx+len*cos(angle);
      *y=cy+len*sin(angle);
      sp=igraph_i_layout_mergegrid_get_sphere(grid, *x, *y, r);
    } while (sp >= 0);

    while (sp < 0 && DIST(*x,*y)<killr) {
      igraph_real_t nx, ny;
      steps++;
      angle=RNG_UNIF(0,2*M_PI);
      len=RNG_UNIF(0, startr/100);
      nx= *x + len * cos(angle);
      ny= *y + len * sin(angle);      
      sp=igraph_i_layout_mergegrid_get_sphere(grid, nx, ny, r);
      if (sp < 0) {
	*x = nx; *y = ny;
      }
    }
  }

/*   fprintf(stderr, "%li ", steps); */
  return 0;
}
