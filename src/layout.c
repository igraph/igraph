/* -*- mode: C -*-  */
/* vim:set ts=2 sw=2 sts=2 et: */
/* 
   IGraph R package.
   Copyright (C) 2003-2014  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA
   
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

#include "igraph_layout.h"
#include "igraph_random.h"
#include "igraph_memory.h"
#include "igraph_iterators.h"
#include "igraph_interface.h"
#include "igraph_adjlist.h"
#include "igraph_progress.h"
#include "igraph_interrupt_internal.h"
#include "igraph_paths.h"
#include "igraph_structural.h"
#include "igraph_visitor.h"
#include "igraph_topology.h"
#include "igraph_components.h"
#include "igraph_types_internal.h"
#include "igraph_dqueue.h"
#include "igraph_arpack.h"
#include "igraph_blas.h"
#include "igraph_centrality.h"
#include "igraph_eigen.h"
#include "config.h"
#include <math.h>
#include "igraph_math.h"


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
 * \param res Pointer to an initialized matrix object. This will
 *        contain the result and will be resized as needed.
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
 * \brief Places the vertices uniformly on a circle, in the order of vertex ids.
 * 
 * \param graph Pointer to an initialized graph object.
 * \param res Pointer to an initialized matrix object. This will
 *        contain the result and will be resized as needed.
 * \param order The order of the vertices on the circle. The vertices
 *        not included here, will be placed at (0,0). Supply
 *        \ref igraph_vss_all() here for all vertices, in the order of
 *        their vertex ids.
 * \return Error code.
 * 
 * Time complexity: O(|V|), the
 * number of vertices. 
 */

int igraph_layout_circle(const igraph_t *graph, igraph_matrix_t *res,
			 igraph_vs_t order) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_integer_t vs_size;
  long int i;
  igraph_vit_t vit;

  IGRAPH_CHECK(igraph_vs_size(graph, &order, &vs_size));

  IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 2));
  igraph_matrix_null(res);

  igraph_vit_create(graph, order, &vit);
  for (i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
    igraph_real_t phi=2*M_PI/vs_size*i;
    int idx = IGRAPH_VIT_GET(vit);
    MATRIX(*res, idx, 0)=cos(phi);
    MATRIX(*res, idx, 1)=sin(phi);
  }
  igraph_vit_destroy(&vit);
  
  return 0;
}

/**
 * \function igraph_layout_star
 * Generate a star-like layout
 * 
 * \param graph The input graph.
 * \param res Pointer to an initialized matrix object. This will
 *        contain the result and will be resized as needed.
 * \param center The id of the vertex to put in the center.
 * \param order A numeric vector giving the order of the vertices 
 *      (including the center vertex!). If a null pointer, then the
 *      vertices are placed in increasing vertex id order.
 * \return Error code.
 * 
 * Time complexity: O(|V|), linear in the number of vertices.
 * 
 * \sa \ref igraph_layout_circle() and other layout generators.
 */

int igraph_layout_star(const igraph_t *graph, igraph_matrix_t *res,
		       igraph_integer_t center, const igraph_vector_t *order) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int c=center;
  long int i;
  igraph_real_t step;
  igraph_real_t phi;

  if (order && igraph_vector_size(order) != no_of_nodes) {
    IGRAPH_ERROR("Invalid order vector length", IGRAPH_EINVAL);
  }
  
  IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 2));

  if (no_of_nodes==1) {
    MATRIX(*res, 0, 0) = MATRIX(*res, 0, 1) = 0.0; 
  } else {
    for (i=0, step=2*M_PI/(no_of_nodes-1), phi=0; 
	 i<no_of_nodes; i++) {
      long int node = order ? (long int) VECTOR(*order)[i] : i;
      if (node != c) { 
	MATRIX(*res, node, 0) = cos(phi);
	MATRIX(*res, node, 1) = sin(phi);
	phi += step;
      } else {
	MATRIX(*res, node, 0) = MATRIX(*res, node, 1) = 0.0;
      }
    }
  }
  
  return 0;
}

/**
 * \function igraph_layout_sphere
 * \brief Places vertices (more or less) uniformly on a sphere.
 *
 * </para><para> 
 * The algorithm was described in the following paper:
 * Distributing many points on a sphere by E.B. Saff and
 * A.B.J. Kuijlaars, \emb Mathematical Intelligencer \eme 19.1 (1997)
 * 5--11.  
 * 
 * \param graph Pointer to an initialized graph object.
 * \param res Pointer to an initialized matrix object. This will
 *        contain the result and will be resized as needed.
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
    IGRAPH_ALLOW_INTERRUPTION();
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
    IGRAPH_ALLOW_INTERRUPTION();
  }
  
  return 0;
}

/**
 * \ingroup layout
 * \function igraph_layout_grid
 * \brief Places the vertices on a regular grid on the plane.
 *
 * \param graph Pointer to an initialized graph object.
 * \param res Pointer to an initialized matrix object. This will
 *        contain the result and will be resized as needed.
 * \param width The number of vertices in a single row of the grid.
 *        When zero or negative, the width of the grid will be the
 *        square root of the number of vertices, rounded up if needed.
 * \return Error code. The current implementation always returns with
 *         success. 
 * 
 * Time complexity: O(|V|), the number of vertices. 
 */
int igraph_layout_grid(const igraph_t *graph, igraph_matrix_t *res, long int width) {
  long int i, no_of_nodes=igraph_vcount(graph);
  igraph_real_t x, y;

  IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 2));  

  if (width <= 0) {
    width = (long int) ceil(sqrt(no_of_nodes));
  }

  x = y = 0;
  for (i = 0; i < no_of_nodes; i++) {
    MATRIX(*res, i, 0) = x++;
    MATRIX(*res, i, 1) = y;
    if (x == width) {
      x = 0; y++;
    }
  }
  
  return 0;
}

/**
 * \ingroup layout
 * \function igraph_layout_grid_3d
 * \brief Places the vertices on a regular grid in the 3D space.
 *
 * \param graph Pointer to an initialized graph object.
 * \param res Pointer to an initialized matrix object. This will
 *        contain the result and will be resized as needed.
 * \param width  The number of vertices in a single row of the grid. When
 *               zero or negative, the width is determined automatically.
 * \param height The number of vertices in a single column of the grid. When
 *               zero or negative, the height is determined automatically.
 *
 * \return Error code. The current implementation always returns with
 *         success. 
 * 
 * Time complexity: O(|V|), the number of vertices. 
 */
int igraph_layout_grid_3d(const igraph_t *graph, igraph_matrix_t *res,
    long int width, long int height) {
  long int i, no_of_nodes=igraph_vcount(graph);
  igraph_real_t x, y, z;

  IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 3));  

  if (width <= 0 && height <= 0) {
    width = height = (long int) ceil(pow(no_of_nodes, 1.0 / 3));
  } else if (width <= 0) {
    width = (long int) ceil(sqrt(no_of_nodes / (double)height));
  } else if (height <= 0) {
    height = (long int) ceil(sqrt(no_of_nodes / (double)width));
  }

  x = y = z = 0;
  for (i = 0; i < no_of_nodes; i++) {
    MATRIX(*res, i, 0) = x++;
    MATRIX(*res, i, 1) = y;
    MATRIX(*res, i, 2) = z;
    if (x == width) {
      x = 0; y++;
      if (y == height) {
        y = 0; z++;
      }
    }
  }
  
  return 0;
}

int igraph_layout_springs(const igraph_t *graph, igraph_matrix_t *res,
			  igraph_real_t mass, igraph_real_t equil, igraph_real_t k,
			  igraph_real_t repeqdis, igraph_real_t kfr, igraph_bool_t repulse) {

  IGRAPH_UNUSED(graph); IGRAPH_UNUSED(res); IGRAPH_UNUSED(mass);
  IGRAPH_UNUSED(equil); IGRAPH_UNUSED(k); IGRAPH_UNUSED(repeqdis);
  IGRAPH_UNUSED(kfr); IGRAPH_UNUSED(repulse);
  IGRAPH_ERROR("Springs layout not implemented", IGRAPH_UNIMPLEMENTED);
  /* TODO */
  return 0;
}

void igraph_i_norm2d(igraph_real_t *x, igraph_real_t *y);

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
 * </para><para>
 * This is a layout generator similar to the Large Graph Layout
 * algorithm and program
 * (http://lgl.sourceforge.net/). But unlike LGL, this
 * version uses a Fruchterman-Reingold style simulated annealing
 * algorithm for placing the vertices. The speedup is achieved by
 * placing the vertices on a grid and calculating the repulsion only
 * for vertices which are closer to each other than a limit. 
 * 
 * \param graph The (initialized) graph object to place.
 * \param res Pointer to an initialized matrix object to hold the
 *   result. It will be resized if needed.
 * \param maxit The maximum number of cooling iterations to perform
 *   for each layout step. A reasonable default is 150.
 * \param maxdelta The maximum length of the move allowed for a vertex
 *   in a single iteration. A reasonable default is the number of
 *   vertices.
 * \param area This parameter gives the area of the square on which
 *   the vertices will be placed. A reasonable default value is the
 *   number of vertices squared.
 * \param coolexp The cooling exponent. A reasonable default value is
 *   1.5.
 * \param repulserad Determines the radius at which vertex-vertex 
 *   repulsion cancels out attraction of adjacent vertices. A
 *   reasonable default value is \p area times the number of vertices.
 * \param cellsize The size of the grid cells, one side of the
 *   square. A reasonable default value is the fourth root of
 *   \p area (or the square root of the number of vertices if \p area
 *   is also left at its default value).
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
  IGRAPH_CHECK(igraph_i_bfs(&mst, (igraph_integer_t) root, IGRAPH_ALL, &vids, 
			    &layers, &parents));
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
  igraph_matrix_scale(res, 1e6);
  
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

    j=(long int) VECTOR(layers)[actlayer];
    for (i=(long int) VECTOR(layers)[actlayer-1]; 
	 i<VECTOR(layers)[actlayer]; i++) {

      long int vid=(long int) VECTOR(vids)[i];
      long int par=(long int) VECTOR(parents)[vid];
      IGRAPH_ALLOW_INTERRUPTION();
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
	igraph_2dgrid_add(&grid, (long int) VECTOR(vids)[j], sx+rx, sy+ry);
	j++; 
      }
    }

    RNG_END();
    
    /*-----------------------------------------*/
    /* Step 2: add the edges of the next layer */
    /*-----------------------------------------*/

    for (j=(long int) VECTOR(layers)[actlayer]; 
	 j<VECTOR(layers)[actlayer+1]; j++) {
      long int vid=(long int) VECTOR(vids)[j];
      long int k;
      IGRAPH_ALLOW_INTERRUPTION();
      IGRAPH_CHECK(igraph_incident(graph, &eids, (igraph_integer_t) vid, 
				   IGRAPH_ALL));
      for (k=0;k<igraph_vector_size(&eids);k++) {
	long int eid=(long int) VECTOR(eids)[k];
	igraph_integer_t from, to;
	igraph_edge(graph, (igraph_integer_t) eid, &from, &to);
	if ((from != vid && igraph_2dgrid_in(&grid, from)) ||
	    (to   != vid && igraph_2dgrid_in(&grid, to))) {
	  igraph_vector_push_back(&edges, eid);
	}
      }
    }

    /*-----------------------------------------*/
    /* Step 3: let the springs spring          */
    /*-----------------------------------------*/
    
    maxchange=epsilon+1;
    while (it < maxit && maxchange > epsilon) {
      long int jj;
      igraph_real_t t=maxdelta*pow((maxit-it)/(double)maxit, coolexp);
      long int vid, nei;

	  IGRAPH_PROGRESS("Large graph layout",
	    100.0*((actlayer-1.0)/(no_of_layers-1.0)+((float)it)/(maxit*(no_of_layers-1.0))),
		0);

      /* init */
      igraph_vector_null(&forcex);
      igraph_vector_null(&forcey);
      maxchange=0;
      
      /* attractive "forces" along the edges */
      for (jj=0; jj<igraph_vector_size(&edges); jj++) {
	igraph_integer_t from, to;
	igraph_real_t xd, yd, dist, force;
	IGRAPH_ALLOW_INTERRUPTION();
	igraph_edge(graph, (igraph_integer_t) VECTOR(edges)[jj], &from, &to);
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
      for (jj=0; jj<VECTOR(layers)[actlayer+1]; jj++) {
	long int vvid=(long int) VECTOR(vids)[jj];
	igraph_real_t fx=VECTOR(forcex)[vvid];
	igraph_real_t fy=VECTOR(forcey)[vvid];
	igraph_real_t ded=sqrt(fx*fx+fy*fy);
	if (ded > t) {
	  ded=t/ded;
	  fx*=ded; fy *=ded;
	}
	igraph_2dgrid_move(&grid, vvid, fx, fy);
	if (fx > maxchange) { maxchange=fx; }
	if (fy > maxchange) { maxchange=fy; }
      }
      it++;
/*       printf("%li iterations, maxchange: %f\n", it, (double)maxchange); */
    }
  }

  IGRAPH_PROGRESS("Large graph layout", 100.0, 0);
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

/* Internal structure for Reingold-Tilford layout */
struct igraph_i_reingold_tilford_vertex {
  long int parent;        /* Parent node index */
  long int level;         /* Level of the node */
  igraph_real_t offset;     /* X offset from parent node */
  long int left_contour;  /* Next left node of the contour
		      of the subtree rooted at this node */
  long int right_contour; /* Next right node of the contour
		      of the subtree rooted at this node */
  igraph_real_t offset_follow_lc;  /* X offset when following the left contour */
  igraph_real_t offset_follow_rc;  /* X offset when following the right contour */
};

int igraph_i_layout_reingold_tilford_postorder(struct igraph_i_reingold_tilford_vertex *vdata,
                                               long int node, long int vcount);
int igraph_i_layout_reingold_tilford_calc_coords(struct igraph_i_reingold_tilford_vertex *vdata,
                                                 igraph_matrix_t *res, long int node,
												 long int vcount, igraph_real_t xpos);

int igraph_i_layout_reingold_tilford(const igraph_t *graph, 
				     igraph_matrix_t *res, 
				     igraph_neimode_t mode, 
				     long int root);
int igraph_i_layout_reingold_tilford(const igraph_t *graph, 
				     igraph_matrix_t *res, igraph_neimode_t mode, 
				     long int root) {
  long int no_of_nodes=igraph_vcount(graph);
  long int i, n, j;
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  igraph_adjlist_t allneis;
  igraph_vector_int_t *neis;
  struct igraph_i_reingold_tilford_vertex *vdata;
    
  IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, 2));
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
  
  IGRAPH_CHECK(igraph_adjlist_init(graph, &allneis, mode));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &allneis);
  
  vdata=igraph_Calloc(no_of_nodes, struct igraph_i_reingold_tilford_vertex);
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
    vdata[i].offset_follow_lc=0.0;
    vdata[i].offset_follow_rc=0.0;
  }
  vdata[root].parent=root;
  vdata[root].level=0;
  MATRIX(*res, root, 1) = 0;
  
  /* Step 1: assign Y coordinates based on BFS and setup parents vector */
  IGRAPH_CHECK(igraph_dqueue_push(&q, root));
  IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
  while (!igraph_dqueue_empty(&q)) {
    long int actnode=(long int) igraph_dqueue_pop(&q);
    long int actdist=(long int) igraph_dqueue_pop(&q);
    neis=igraph_adjlist_get(&allneis, actnode);
    n=igraph_vector_int_size(neis);
    for (j=0; j<n; j++) {
      long int neighbor=(long int) VECTOR(*neis)[j];
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
  igraph_adjlist_destroy(&allneis);
  igraph_free(vdata);
  IGRAPH_FINALLY_CLEAN(3);
  
  IGRAPH_PROGRESS("Reingold-Tilford tree layout", 100.0, NULL);
  
  return 0;
}

int igraph_i_layout_reingold_tilford_calc_coords(struct igraph_i_reingold_tilford_vertex *vdata,
                                                 igraph_matrix_t *res, long int node,
												 long int vcount, igraph_real_t xpos) {
  long int i;
  MATRIX(*res, node, 0) = xpos;
  for (i=0; i<vcount; i++) {
    if (i == node) continue;
    if (vdata[i].parent == node) {
      igraph_i_layout_reingold_tilford_calc_coords(vdata, res, i, vcount,
						   xpos+vdata[i].offset);
    }
  }
  return 0;
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
  /*printf("Visited node %d and arranged its subtrees\n", node);*/
  for (i=0, j=0; i<vcount; i++) {
    if (i == node) continue;
    if (vdata[i].parent == node) {
      /*printf("  Placing child %d on level %d\n", i, vdata[i].level);*/
      if (leftroot >= 0) {
        /* Now we will follow the right contour of leftroot and the
         * left contour of the subtree rooted at i */
        long lnode, rnode;
        igraph_real_t loffset, roffset, minsep, rootsep;
        lnode = leftroot; rnode = i;
        minsep = 1;
        rootsep = vdata[leftroot].offset + minsep;
        loffset = 0; roffset = minsep;
        /*printf("    Contour: [%d, %d], offsets: [%lf, %lf], rootsep: %lf\n",
               lnode, rnode, loffset, roffset, rootsep);*/
        while ((lnode >= 0) && (rnode >= 0)) {
          /* Step to the next level on the right contour of the left subtree */
          if (vdata[lnode].right_contour >= 0) {
            loffset += vdata[lnode].offset_follow_rc;
            lnode = vdata[lnode].right_contour;
          } else {
            /* Left subtree ended there. The right contour of the left subtree
             * will continue to the next step on the right subtree. */
            if (vdata[rnode].left_contour >= 0) {
              /*printf("      Left subtree ended, continuing left subtree's left and right contour on right subtree (node %ld)\n", vdata[rnode].left_contour);*/
              vdata[lnode].left_contour = vdata[rnode].left_contour;
              vdata[lnode].right_contour = vdata[rnode].left_contour;
              vdata[lnode].offset_follow_lc = vdata[lnode].offset_follow_rc =
                (roffset-loffset)+vdata[rnode].offset_follow_lc;
              /*printf("      vdata[lnode].offset_follow_* = %.4f\n", vdata[lnode].offset_follow_lc);*/
            }
            lnode = -1;
          }
          /* Step to the next level on the left contour of the right subtree */
          if (vdata[rnode].left_contour >= 0) {
            roffset += vdata[rnode].offset_follow_lc;
            rnode = vdata[rnode].left_contour;
          } else {
            /* Right subtree ended here. The left contour of the right
             * subtree will continue to the next step on the left subtree.
             * Note that lnode has already been advanced here */
            if (lnode >= 0) {
              /*printf("      Right subtree ended, continuing right subtree's left and right contour on left subtree (node %ld)\n", lnode);*/
              vdata[rnode].left_contour = lnode;
              vdata[rnode].right_contour = lnode;
              vdata[rnode].offset_follow_lc = vdata[rnode].offset_follow_rc =
                (loffset-roffset);  /* loffset has also been increased earlier */
              /*printf("      vdata[rnode].offset_follow_* = %.4f\n", vdata[rnode].offset_follow_lc);*/
            }
            rnode = -1;
          }
          /*printf("    Contour: [%d, %d], offsets: [%lf, %lf], rootsep: %lf\n", 
                 lnode, rnode, loffset, roffset, rootsep);*/
      
          /* Push subtrees away if necessary */
          if ((lnode >= 0) && (rnode >= 0) && (roffset - loffset < minsep)) {
            /*printf("    Pushing right subtree away by %lf\n", minsep-roffset+loffset);*/
            rootsep += minsep-roffset+loffset;
            roffset = loffset+minsep;
          }
        }

        /*printf("  Offset of subtree with root node %d will be %lf\n", i, rootsep);*/
        vdata[i].offset = rootsep;
        vdata[node].right_contour = i;
        vdata[node].offset_follow_rc = rootsep;
        avg = (avg*j)/(j+1) + rootsep/(j+1);
        leftrootidx=j;
        leftroot=i;
      } else {
        leftrootidx=j;
        leftroot=i;
        vdata[node].left_contour=i;
        vdata[node].right_contour=i;
        vdata[node].offset_follow_lc = 0.0;
        vdata[node].offset_follow_rc = 0.0;
        avg = vdata[i].offset; 
      }
      j++;
    }
  }
  /*printf("Shifting node to be centered above children. Shift amount: %lf\n", avg);*/
  vdata[node].offset_follow_lc -= avg;
  vdata[node].offset_follow_rc -= avg;
  for (i=0, j=0; i<vcount; i++) {
    if (i == node) continue;
    if (vdata[i].parent == node) vdata[i].offset -= avg;
  }
  
  return 0;
}

/**
 * \function igraph_layout_reingold_tilford
 * \brief Reingold-Tilford layout for tree graphs
 * 
 * </para><para>
 * Arranges the nodes in a tree where the given node is used as the root.
 * The tree is directed downwards and the parents are centered above its
 * children. For the exact algorithm, see:
 * 
 * </para><para>
 * Reingold, E and Tilford, J: Tidier drawing of trees.
 * IEEE Trans. Softw. Eng., SE-7(2):223--228, 1981
 *
 * </para><para>
 * If the given graph is not a tree, a breadth-first search is executed
 * first to obtain a possible spanning tree.
 * 
 * \param graph The graph object. 
 * \param res The result, the coordinates in a matrix. The parameter
 *   should point to an initialized matrix object and will be resized.
 * \param mode Specifies which edges to consider when building the tree.
 *   If it is \c IGRAPH_OUT then only the outgoing, if it is \c IGRAPH_IN
 *   then only the incoming edges of a parent are considered. If it is 
 *   \c IGRAPH_ALL then all edges are used (this was the behavior in 
 *   igraph 0.5 and before). This parameter also influences how the root 
 *   vertices are calculated, if they are not given. See the \p roots parameter.
 * \param roots The index of the root vertex or root vertices. 
 *   If this is a non-empty vector then the supplied vertex ids are used 
 *   as the roots of the trees (or a single tree if the graph is connected).
 *   If it is a null pointer of a pointer to an empty vector, then the root
 *   vertices are automatically calculated based on topological sorting,
 *   performed with the opposite mode than the \p mode argument.
 *   After the vertices have been sorted, one is selected from each component.
 * \param rootlevel This argument can be useful when drawing forests which are 
 *   not trees (i.e. they are unconnected and have tree components). It specifies 
 *   the level of the root vertices for every tree in the forest. It is only
 *   considered if not a null pointer and the \p roots argument is also given 
 *   (and it is not a null pointer of an empty vector).
 * \return Error code.
 *
 * Added in version 0.2.
 * 
 * \sa \ref igraph_layout_reingold_tilford_circular().
 * 
 * \example examples/simple/igraph_layout_reingold_tilford.c
 */

int igraph_layout_reingold_tilford(const igraph_t *graph,
				   igraph_matrix_t *res, 
				   igraph_neimode_t mode,
				   const igraph_vector_t *roots,
				   const igraph_vector_t *rootlevel) {

  long int no_of_nodes_orig=igraph_vcount(graph);
  long int no_of_nodes=no_of_nodes_orig;
  long int real_root;
  igraph_t extended;
  const igraph_t *pextended=graph;
  igraph_vector_t myroots;
  const igraph_vector_t *proots=roots;
  igraph_neimode_t mode2;
  
  /* TODO: possible speedup could be achieved if we use a table for storing
   * the children of each node in the tree. (Now the implementation uses a
   * single array containing the parent of each node and a node's children
   * are determined by looking for other nodes that have this node as parent)
   */

  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  if (mode==IGRAPH_IN) {
    mode2=IGRAPH_OUT;
  } else if (mode==IGRAPH_OUT) {
    mode2=IGRAPH_IN;
  } else {
    mode2=mode;
  }

  if ( (!roots || igraph_vector_size(roots)==0) && 
       rootlevel && igraph_vector_size(rootlevel) != 0 ) {
    IGRAPH_WARNING("Reingold-Tilford layout: 'rootlevel' ignored");
  }

  /* ----------------------------------------------------------------------- */
  /* If root vertices are not given, then do a topological sort and take 
     the last element from every component for directed graphs, or select the
     vertex with the maximum degree from each component for undirected graphs */

  if (!roots || igraph_vector_size(roots)==0) {
    
    igraph_vector_t order, membership;
    igraph_integer_t no_comps;
    long int i, noseen=0;

    IGRAPH_VECTOR_INIT_FINALLY(&myroots, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&order, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&membership, no_of_nodes);

    if (igraph_is_directed(graph) && mode != IGRAPH_ALL) {
      IGRAPH_CHECK(igraph_topological_sorting(graph, &order, mode2));
      IGRAPH_CHECK(igraph_clusters(graph, &membership, /*csize=*/ 0, 
        &no_comps, IGRAPH_WEAK));
    } else {
      IGRAPH_CHECK(igraph_sort_vertex_ids_by_degree(graph, &order,
            igraph_vss_all(), IGRAPH_ALL, 0, IGRAPH_ASCENDING, 0));
      IGRAPH_CHECK(igraph_clusters(graph, &membership, /*csize=*/ 0, 
        &no_comps, IGRAPH_WEAK));
    }
    
    IGRAPH_CHECK(igraph_vector_resize(&myroots, no_comps));
    igraph_vector_null(&myroots);
    proots=&myroots;
    for (i=no_of_nodes-1; noseen < no_comps && i>=0; i--) {
      long int v=(long int) VECTOR(order)[i];
      long int mem=(long int) VECTOR(membership)[v];
      if (VECTOR(myroots)[mem]==0) {
        noseen += 1;
        VECTOR(myroots)[mem]=v+1;
      }
    }

    for (i=0; i<no_comps; i++) {
      VECTOR(myroots)[i] -= 1;
    }

    igraph_vector_destroy(&membership);
    igraph_vector_destroy(&order);
    IGRAPH_FINALLY_CLEAN(2);

  } else if (rootlevel && igraph_vector_size(rootlevel) > 0 &&
	     igraph_vector_size(roots) > 1) {

  /* ----------------------------------------------------------------------- */
  /* Many roots were given to us, check 'rootlevel' */

    long int plus_levels=0;
    long int i;
    
    if (igraph_vector_size(roots) != igraph_vector_size(rootlevel)) {
      IGRAPH_ERROR("Reingold-Tilford: 'roots' and 'rootlevel' lengths differ",
		   IGRAPH_EINVAL);
    }

    /* check if there is one which is not zero */
    for (i=0; i<igraph_vector_size(roots); i++) {
      plus_levels += VECTOR(*rootlevel)[i];
    }

    /* make copy of graph, add vertices/edges */
    if (plus_levels != 0) {
      igraph_vector_t newedges;
      long int edgeptr=0;

      pextended=&extended;
      IGRAPH_CHECK(igraph_copy(&extended, graph));
      IGRAPH_FINALLY(igraph_destroy, &extended);
      IGRAPH_CHECK(igraph_add_vertices(&extended, 
				       (igraph_integer_t) plus_levels, 0));

      IGRAPH_VECTOR_INIT_FINALLY(&newedges, plus_levels*2);

      for (i=0; i<igraph_vector_size(roots); i++) {
	long int rl=(long int) VECTOR(*rootlevel)[i];
	long int rn=(long int) VECTOR(*roots)[i];
	long int j;

	if (rl==0) { continue; }

	VECTOR(newedges)[edgeptr++] = no_of_nodes;
	VECTOR(newedges)[edgeptr++] = rn;
	
	for (j=0; j<rl-1; j++) {
	  VECTOR(newedges)[edgeptr++] = no_of_nodes+1;
	  VECTOR(newedges)[edgeptr++] = no_of_nodes;
	  no_of_nodes++;
	}

	VECTOR(*roots)[i]=no_of_nodes++;
      }
      
      /* Opposite direction if mode=="in" */
      if (mode==IGRAPH_IN) {
	for (i=0; i<igraph_vector_size(&newedges); i+=2) {
	  igraph_real_t tmp=VECTOR(newedges)[i];
	  VECTOR(newedges)[i] = VECTOR(newedges)[i+1];
	  VECTOR(newedges)[i+1] = tmp;
	}
      }

      IGRAPH_CHECK(igraph_add_edges(&extended, &newedges, 0));
      igraph_vector_destroy(&newedges);
      IGRAPH_FINALLY_CLEAN(1);
    }
  }
  
  /* ----------------------------------------------------------------------- */
  /* Ok, we have the root vertices now. 
     Possibly we need to add some edges to make the graph connected. */

  if (igraph_vector_size(proots)==1) {
    real_root=(long int) VECTOR(*proots)[0];
    if (real_root<0 || real_root>=no_of_nodes) {
      IGRAPH_ERROR("invalid vertex id", IGRAPH_EINVVID);
    }
  } else {
    igraph_vector_t newedges;
    long int no_of_newedges=igraph_vector_size(proots);
    long int i;
    real_root=no_of_nodes;
    
    /* Make copy if needed */
    if (pextended == graph) {
      pextended=&extended;
      IGRAPH_CHECK(igraph_copy(&extended, graph));
      IGRAPH_FINALLY(igraph_destroy, &extended);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&newedges, no_of_newedges*2);
    IGRAPH_CHECK(igraph_add_vertices(&extended, 1, 0));
    for (i=0; i<no_of_newedges; i++) {
      VECTOR(newedges)[2*i] = no_of_nodes;
      VECTOR(newedges)[2*i+1] = VECTOR(*proots)[i];
    }
    IGRAPH_CHECK(igraph_add_edges(&extended, &newedges, 0));
    
    no_of_nodes++;
    igraph_vector_destroy(&newedges);
    IGRAPH_FINALLY_CLEAN(1);
  }

  /* ----------------------------------------------------------------------- */
  /* Layout */

  IGRAPH_CHECK(igraph_i_layout_reingold_tilford(pextended, res, mode, real_root));
  
  /* Remove the new vertices from the layout */
  if (no_of_nodes != no_of_nodes_orig) {
    if (no_of_nodes-1 == no_of_nodes_orig) {
      IGRAPH_CHECK(igraph_matrix_remove_row(res, no_of_nodes_orig));
    } else {
      igraph_matrix_t tmp;
      long int i;
      IGRAPH_MATRIX_INIT_FINALLY(&tmp, no_of_nodes_orig, 2);
      for (i=0; i<no_of_nodes_orig; i++) {
	MATRIX(tmp, i, 0) = MATRIX(*res, i, 0);
	MATRIX(tmp, i, 1) = MATRIX(*res, i, 1);
      }
      IGRAPH_CHECK(igraph_matrix_update(res, &tmp));
      igraph_matrix_destroy(&tmp);
      IGRAPH_FINALLY_CLEAN(1);
    }
  }
  
  if (pextended != graph) {
    igraph_destroy(&extended);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  /* Remove the roots vector if it was created by us */
  if (proots != roots) {
    igraph_vector_destroy(&myroots);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  return 0;
}

/** 
 * \function igraph_layout_reingold_tilford_circular
 * \brief Circular Reingold-Tilford layout for trees
 * 
 * </para><para>
 * This layout is almost the same as \ref igraph_layout_reingold_tilford(), but 
 * the tree is drawn in a circular way, with the root vertex in the center.
 * 
 * \param graph The graph object.
 * \param res The result, the coordinates in a matrix. The parameter
 *   should point to an initialized matrix object and will be resized.
 * \param mode Specifies which edges to consider when building the tree.
 *   If it is \c IGRAPH_OUT then only the outgoing, if it is \c IGRAPH_IN
 *   then only the incoming edges of a parent are considered. If it is 
 *   \c IGRAPH_ALL then all edges are used (this was the behavior in 
 *   igraph 0.5 and before). This parameter also influences how the root 
 *   vertices are calculated, if they are not given. See the \p roots parameter.
 * \param roots The index of the root vertex or root vertices. 
 *   If this is a non-empty vector then the supplied vertex ids are used 
 *   as the roots of the trees (or a single tree if the graph is connected).
 *   If it is a null pointer of a pointer to an empty vector, then the root
 *   vertices are automatically calculated based on topological sorting,
 *   performed with the opposite mode than the \p mode argument.
 *   After the vertices have been sorted, one is selected from each component.
 * \param rootlevel This argument can be useful when drawing forests which are 
 *   not trees (i.e. they are unconnected and have tree components). It specifies 
 *   the level of the root vertices for every tree in the forest. It is only
 *   considered if not a null pointer and the \p roots argument is also given 
 *   (and it is not a null pointer of an empty vector). Note that if you supply
 *   a null pointer here and the graph has multiple components, all of the root
 *   vertices will be mapped to the origin of the coordinate system, which does
 *   not really make sense.
 * \return Error code.
 *
 * \sa \ref igraph_layout_reingold_tilford().
 */

int igraph_layout_reingold_tilford_circular(const igraph_t *graph,
					    igraph_matrix_t *res,
					    igraph_neimode_t mode,
					    const igraph_vector_t *roots,
					    const igraph_vector_t *rootlevel) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int i;
  igraph_real_t ratio=2*M_PI*(no_of_nodes-1.0)/no_of_nodes;
  igraph_real_t minx, maxx;

  IGRAPH_CHECK(igraph_layout_reingold_tilford(graph, res, mode, roots, rootlevel));

  if (no_of_nodes == 0) return 0;

  minx = maxx = MATRIX(*res, 0, 0);
  for (i=1; i<no_of_nodes; i++) {
    if (MATRIX(*res, i, 0) > maxx) maxx=MATRIX(*res, i, 0);
	if (MATRIX(*res, i, 0) < minx) minx=MATRIX(*res, i, 0);
  }
  if (maxx > minx)
    ratio /= (maxx-minx);
  for (i=0; i<no_of_nodes; i++) {
    igraph_real_t phi=(MATRIX(*res, i, 0)-minx)*ratio;
    igraph_real_t r=MATRIX(*res, i, 1);
    MATRIX(*res, i, 0) = r*cos(phi);
    MATRIX(*res, i, 1) = r*sin(phi);
  }
  
  return 0;
}

#define COULOMBS_CONSTANT 8987500000.0


igraph_real_t igraph_i_distance_between(const igraph_matrix_t *c, long int a, 
					long int b);

int igraph_i_determine_electric_axal_forces(const igraph_matrix_t *pos,
					    igraph_real_t *x,
					    igraph_real_t *y,
					    igraph_real_t directed_force,
					    igraph_real_t distance,
					    long int other_node,
					    long int this_node);

int igraph_i_apply_electrical_force(const igraph_matrix_t *pos,
				    igraph_vector_t *pending_forces_x,
				    igraph_vector_t *pending_forces_y,
				    long int other_node, long int this_node,
				    igraph_real_t node_charge, 
				    igraph_real_t distance);

int igraph_i_determine_spring_axal_forces(const igraph_matrix_t *pos,
					  igraph_real_t *x, igraph_real_t *y,
					  igraph_real_t directed_force,
					  igraph_real_t distance,
					  int spring_length,
					  long int other_node, 
					  long int this_node);

int igraph_i_apply_spring_force(const igraph_matrix_t *pos, 
				igraph_vector_t *pending_forces_x,
				igraph_vector_t *pending_forces_y,
				long int other_node,
				long int this_node, int spring_length,
				igraph_real_t spring_constant);

int igraph_i_move_nodes(igraph_matrix_t *pos, 
			const igraph_vector_t *pending_forces_x, 
			const igraph_vector_t *pending_forces_y,
			igraph_real_t node_mass,
			igraph_real_t max_sa_movement);

igraph_real_t igraph_i_distance_between(const igraph_matrix_t *c, long int a, 
					long int b) {
  igraph_real_t diffx=MATRIX(*c, a, 0)-MATRIX(*c, b, 0);
  igraph_real_t diffy=MATRIX(*c, a, 1)-MATRIX(*c, b, 1);
  return sqrt( diffx*diffx + diffy*diffy );
} 

int igraph_i_determine_electric_axal_forces(const igraph_matrix_t *pos,
					    igraph_real_t *x,
					    igraph_real_t *y,
					    igraph_real_t directed_force,
					    igraph_real_t distance,
					    long int other_node,
					    long int this_node) {
  
  // We know what the directed force is.  We now need to translate it
  // into the appropriate x and y components.
  // First, assume: 
  //                 other_node
  //                    /|
  //  directed_force  /  |
  //                /    | y
  //              /______|
  //    this_node     x         
  //
  // other_node.x > this_node.x
  // other_node.y > this_node.y
  // the force will be on this_node away from other_node
  
  // the proportion (distance/y_distance) is equal to the proportion
  // (directed_force/y_force), as the two triangles are similar.
  // therefore, the magnitude of y_force = (directed_force*y_distance)/distance
  // the sign of y_force is negative, away from other_node
  
  igraph_real_t x_distance, y_distance;
  y_distance = MATRIX(*pos, other_node, 1)-MATRIX(*pos, this_node, 1);
  if (y_distance < 0) { y_distance = -y_distance; }
  *y = -1 * ((directed_force * y_distance) / distance);

  // the x component works in exactly the same way.
  x_distance = MATRIX(*pos, other_node, 0)-MATRIX(*pos, this_node, 0);
  if (x_distance < 0) { x_distance = -x_distance; }
  *x = -1 * ((directed_force * x_distance) / distance);
  
   // Now we need to reverse the polarity of our answers based on the falsness
   // of our assumptions.
  if (MATRIX(*pos, other_node, 0) < MATRIX(*pos, this_node, 0)) {
    *x = *x * -1;
  }
  if (MATRIX(*pos, other_node, 1) < MATRIX(*pos, this_node, 1)) {
    *y = *y * -1;
  }
  
  return 0;
}
  
int igraph_i_apply_electrical_force(const igraph_matrix_t *pos,
				    igraph_vector_t *pending_forces_x,
				    igraph_vector_t *pending_forces_y,
				    long int other_node, long int this_node,
				    igraph_real_t node_charge, 
				    igraph_real_t distance) {

  igraph_real_t directed_force = COULOMBS_CONSTANT * 
    ((node_charge * node_charge)/(distance * distance));
  
  igraph_real_t x_force, y_force;
  igraph_i_determine_electric_axal_forces(pos, &x_force, &y_force, 
					  directed_force, distance, 
					  other_node, this_node);

  VECTOR(*pending_forces_x)[this_node] += x_force;
  VECTOR(*pending_forces_y)[this_node] += y_force;
  VECTOR(*pending_forces_x)[other_node] -= x_force;
  VECTOR(*pending_forces_y)[other_node] -= y_force;

  return 0;
}

int igraph_i_determine_spring_axal_forces(const igraph_matrix_t *pos,
					  igraph_real_t *x, igraph_real_t *y,
					  igraph_real_t directed_force,
					  igraph_real_t distance,
					  int spring_length,
					  long int other_node, long int this_node) {

  // if the spring is just the right size, the forces will be 0, so we can
  // skip the computation.
  //
  // if the spring is too long, our forces will be identical to those computed
  // by determine_electrical_axal_forces() (this_node will be pulled toward
  // other_node).
  //
  // if the spring is too short, our forces will be the opposite of those
  // computed by determine_electrical_axal_forces() (this_node will be pushed
  // away from other_node)
  //
  // finally, since both nodes are movable, only one-half of the total force
  // should be applied to each node, so half the forces for our answer.
  
  if (distance == spring_length) {
    *x = 0.0;
    *y = 0.0;
  } else {
    igraph_i_determine_electric_axal_forces(pos, x, y, directed_force, distance, 
					    other_node, this_node);
    if (distance < spring_length) {
      *x = -1 * *x;
      *y = -1 * *y;
    }
    *x = 0.5 * *x;
    *y = 0.5 * *y;
  }
  
  return 0;
}

int igraph_i_apply_spring_force(const igraph_matrix_t *pos, 
				igraph_vector_t *pending_forces_x,
				igraph_vector_t *pending_forces_y,
				long int other_node,
				long int this_node, int spring_length,
				igraph_real_t spring_constant) {

  // determined using Hooke's Law:
  //   force = -kx
  // where:
  //   k = spring constant
  //   x = displacement from ideal length in meters
  
  igraph_real_t distance, displacement, directed_force, x_force, y_force;
  distance = igraph_i_distance_between(pos, other_node, this_node);
  // let's protect ourselves from division by zero by ignoring two nodes that
  // happen to be in the same place.  Since we separate all nodes before we
  // work on any of them, this will only happen in extremely rare circumstances,
  // and when it does, electrical force will probably push one or both of them
   // one way or another anyway.
  if (distance == 0.0) {
    return 0;
  }
  
  displacement = distance - spring_length;
  if (displacement < 0) {
    displacement = -displacement;
  }
  directed_force = -1 * spring_constant * displacement;
  // remember, this is force directed away from the spring;
  // a negative number is back towards the spring (or, in our case, back towards
   // the other node)
  
  // get the force that should be applied to >this< node
  igraph_i_determine_spring_axal_forces(pos, &x_force, &y_force, 
					directed_force, distance, spring_length,
					other_node, this_node);
  
  VECTOR(*pending_forces_x)[this_node] += x_force;
  VECTOR(*pending_forces_y)[this_node] += y_force;
  VECTOR(*pending_forces_x)[other_node] -= x_force;
  VECTOR(*pending_forces_y)[other_node] -= y_force;

  return 0;
}

int igraph_i_move_nodes(igraph_matrix_t *pos, 
			const igraph_vector_t *pending_forces_x, 
			const igraph_vector_t *pending_forces_y,
			igraph_real_t node_mass,
			igraph_real_t max_sa_movement) {
  
  // Since each iteration is isolated, time is constant at 1.
  // Therefore:
  //   Force effects acceleration.
  //   acceleration (d(velocity)/time) = velocity
  //   velocity (d(displacement)/time) = displacement
  //   displacement = acceleration
  
  // determined using Newton's second law:
  //   sum(F) = ma
  // therefore:
  //   acceleration = force / mass
  //   velocity     = force / mass
  //   displacement = force / mass

  long int this_node, no_of_nodes=igraph_vector_size(pending_forces_x);

  for (this_node=0; this_node < no_of_nodes; this_node++) {

    igraph_real_t x_movement, y_movement;

    x_movement = VECTOR(*pending_forces_x)[this_node] / node_mass;
    if (x_movement > max_sa_movement) {
      x_movement = max_sa_movement;
    } else if (x_movement < -max_sa_movement) {
      x_movement = -max_sa_movement;
    }

    y_movement = VECTOR(*pending_forces_y)[this_node] / node_mass;
    if (y_movement > max_sa_movement) {
      y_movement = max_sa_movement;
    } else if (y_movement < -max_sa_movement) {
      y_movement = -max_sa_movement;
    }

    MATRIX(*pos, this_node, 0) += x_movement;
    MATRIX(*pos, this_node, 1) += y_movement;

  }
  return 0;
}

/**
 * \function igraph_layout_graphopt
 * \brief Optimizes vertex layout via the graphopt algorithm.
 * 
 * </para><para>
 * This is a port of the graphopt layout algorithm by Michael Schmuhl.
 * graphopt version 0.4.1 was rewritten in C and the support for 
 * layers was removed (might be added later) and a code was a bit 
 * reorganized to avoid some unnecessary steps is the node charge (see below) 
 * is zero.
 * 
 * </para><para>
 * graphopt uses physical analogies for defining attracting and repelling 
 * forces among the vertices and then the physical system is simulated 
 * until it reaches an equilibrium. (There is no simulated annealing or 
 * anything like that, so a stable fixed point is not guaranteed.)
 * 
 * </para><para>
 * See also http://www.schmuhl.org/graphopt/ for the original graphopt.
 * \param graph The input graph.
 * \param res Pointer to an initialized matrix, the result will be stored here
 *    and its initial contents is used the starting point of the simulation
 *    if the \p use_seed argument is true. Note that in this case the 
 *    matrix should have the proper size, otherwise a warning is issued and 
 *    the supplied values are ignored. If no starting positions are given 
 *    (or they are invalid) then a random staring position is used. 
 *    The matrix will be resized if needed.
 * \param niter Integer constant, the number of iterations to perform.
 *    Should be a couple of hundred in general. If you have a large graph 
 *    then you might want to only do a few iterations and then check the 
 *    result. If it is not good enough you can feed it in again in 
 *    the \p res argument. The original graphopt default if 500.
 * \param node_charge The charge of the vertices, used to calculate electric
 *    repulsion. The original graphopt default is 0.001.
 * \param node_mass The mass of the vertices, used for the spring forces.
 *    The original graphopt defaults to 30.
 * \param spring_length The length of the springs, an integer number.
 *    The original graphopt defaults to zero.
 * \param spring_constant The spring constant, the original graphopt defaults 
 *    to one.
 * \param max_sa_movement Real constant, it gives the maximum amount of movement 
 *    allowed in a single step along a single axis. The original graphopt 
 *    default is 5.
 * \param use_seed Logical scalar, whether to use the positions in \p res as
 *    a starting configuration. See also \p res above.
 * \return Error code.
 * 
 * Time complexity: O(n (|V|^2+|E|) ), n is the number of iterations, 
 * |V| is the number of vertices, |E| the number
 * of edges. If \p node_charge is zero then it is only O(n|E|).
 */

int igraph_layout_graphopt(const igraph_t *graph, igraph_matrix_t *res, 
			   igraph_integer_t niter,
			   igraph_real_t node_charge, igraph_real_t node_mass,
			   igraph_real_t spring_length,
			   igraph_real_t spring_constant, 
			   igraph_real_t max_sa_movement,
			   igraph_bool_t use_seed) {

  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  int my_spring_length=(int) spring_length;
  igraph_vector_t pending_forces_x, pending_forces_y;
  /* Set a flag to calculate (or not) the electrical forces that the nodes */
  /* apply on each other based on if both node types' charges are zero. */
  igraph_bool_t apply_electric_charges= (node_charge!=0);
  
  long int this_node, other_node, edge;
  igraph_real_t distance;
  long int i;

  IGRAPH_VECTOR_INIT_FINALLY(&pending_forces_x, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&pending_forces_y, no_of_nodes);
  
  if (use_seed) {
    if (igraph_matrix_nrow(res) != no_of_nodes ||
	igraph_matrix_ncol(res) != 2) {
      IGRAPH_WARNING("Invalid size for initial matrix, starting from random layout");
      IGRAPH_CHECK(igraph_layout_random(graph, res));
    }
  } else {
    IGRAPH_CHECK(igraph_layout_random(graph, res));
  }

  IGRAPH_PROGRESS("Graphopt layout", 0, NULL);
  for(i=niter;i>0;i--) {
    /* Report progress in approx. every 100th step */
    if (i%10 == 0) {
      IGRAPH_PROGRESS("Graphopt layout", 100.0-100.0*i/niter, NULL);
    }
    
    /* Clear pending forces on all nodes */
    igraph_vector_null(&pending_forces_x);
    igraph_vector_null(&pending_forces_y);
    
    // Apply electrical force applied by all other nodes
    if (apply_electric_charges) {
      // Iterate through all nodes
      for (this_node = 0; this_node < no_of_nodes; this_node++) {
	IGRAPH_ALLOW_INTERRUPTION();
	for (other_node = this_node + 1; 
	     other_node < no_of_nodes; 
	     other_node++) {
	  distance = igraph_i_distance_between(res, this_node, other_node);
	  // let's protect ourselves from division by zero by ignoring
	  // two nodes that happen to be in the same place.  Since we 
	  // separate all nodes before we work on any of them, this 
	  // will only happen in extremely rare circumstances, and when
	  // it does, springs will probably pull them apart anyway.
	  // also, if we are more than 50 away, the electric force 
	  // will be negligible.  
	  // ***** may not always be desirable ****
	  if ((distance != 0.0) && (distance < 500.0)) {
	    //	  if (distance != 0.0) {
	    // Apply electrical force from node(counter2) on 
	    // node(counter)
	    igraph_i_apply_electrical_force(res, &pending_forces_x, 
					    &pending_forces_y, 
					    other_node, this_node, 
					    node_charge,
					    distance);
	  }
	}
      }
    }
      
    // Apply force from springs
    for (edge = 0; edge < no_of_edges; edge++) {
      long int tthis_node=IGRAPH_FROM(graph, edge);
      long int oother_node=IGRAPH_TO(graph, edge);
      // Apply spring force on both nodes
      igraph_i_apply_spring_force(res, &pending_forces_x, &pending_forces_y,
				  oother_node, tthis_node, my_spring_length, 
				  spring_constant);
    }
  
    // Effect the movement of the nodes based on all pending forces
    igraph_i_move_nodes(res, &pending_forces_x, &pending_forces_y, node_mass,
			max_sa_movement);
  }
  IGRAPH_PROGRESS("Graphopt layout", 100, NULL);

  igraph_vector_destroy(&pending_forces_y);
  igraph_vector_destroy(&pending_forces_x);
  IGRAPH_FINALLY_CLEAN(2);
  
  return 0;
}

int igraph_i_layout_merge_dla(igraph_i_layout_mergegrid_t *grid, 
			      long int actg, igraph_real_t *x, igraph_real_t *y, igraph_real_t r,
			      igraph_real_t cx, igraph_real_t cy, igraph_real_t startr, 
			      igraph_real_t killr);

int igraph_i_layout_sphere_2d(igraph_matrix_t *coords, igraph_real_t *x, 
			      igraph_real_t *y, igraph_real_t *r);
int igraph_i_layout_sphere_3d(igraph_matrix_t *coords, igraph_real_t *x,
			      igraph_real_t *y, igraph_real_t *z,
			      igraph_real_t *r);

/**
 * \function igraph_layout_merge_dla
 * \brief Merge multiple layouts by using a DLA algorithm
 * 
 * </para><para>
 * First each layout is covered by a circle. Then the layout of the
 * largest graph is placed at the origin. Then the other layouts are
 * placed by the DLA algorithm, larger ones first and smaller ones
 * last.
 * \param thegraphs Pointer vector containing the graph object of
 *        which the layouts will be merged.
 * \param coords Pointer vector containing matrix objects with the 2d
 *        layouts of the graphs in \p thegraphs.
 * \param res Pointer to an initialized matrix object, the result will
 *        be stored here. It will be resized if needed.
 * \return Error code.
 * 
 * Added in version 0.2. This function is experimental.
 * 
 * </para><para>
 * Time complexity: TODO.
 */

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

  /* Graphs are currently not used, only the coordinates */
  IGRAPH_UNUSED(thegraphs);
  
  IGRAPH_VECTOR_INIT_FINALLY(&sizes, graphs);
  IGRAPH_VECTOR_INIT_FINALLY(&x, graphs);
  IGRAPH_VECTOR_INIT_FINALLY(&y, graphs);
  IGRAPH_VECTOR_INIT_FINALLY(&r, graphs);
  IGRAPH_VECTOR_INIT_FINALLY(&nx, graphs);
  IGRAPH_VECTOR_INIT_FINALLY(&ny, graphs);
  IGRAPH_VECTOR_INIT_FINALLY(&nr, graphs);

  RNG_BEGIN();
  
  for (i=0; i<igraph_vector_ptr_size(coords); i++) {
    igraph_matrix_t *mat=VECTOR(*coords)[i];
    long int size=igraph_matrix_nrow(mat);

    if (igraph_matrix_ncol(mat) != 2) {
      IGRAPH_ERROR("igraph_layout_merge_dla works for 2D layouts only",
                   IGRAPH_EINVAL);
    }

    IGRAPH_ALLOW_INTERRUPTION();
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
  actg=(long int) VECTOR(sizes)[jpos++];
  igraph_i_layout_merge_place_sphere(&grid, 0, 0, VECTOR(r)[actg], actg);
  
  IGRAPH_PROGRESS("Merging layouts via DLA", 0.0, NULL);
  while (jpos<graphs) {
    IGRAPH_ALLOW_INTERRUPTION();
/*     fprintf(stderr, "comp: %li", jpos); */
    IGRAPH_PROGRESS("Merging layouts via DLA", (100.0*jpos)/graphs, NULL);
    
    actg=(long int) VECTOR(sizes)[jpos++];
    /* 2. random walk, TODO: tune parameters */
    igraph_i_layout_merge_dla(&grid, actg, 
			      igraph_vector_e_ptr(&x, actg),
			      igraph_vector_e_ptr(&y, actg), 
			      VECTOR(r)[actg], 0, 0,
			      maxx, maxx+5);
    
    /* 3. place sphere */
    igraph_i_layout_merge_place_sphere(&grid, VECTOR(x)[actg], VECTOR(y)[actg],
				       VECTOR(r)[actg], actg);
  }
  IGRAPH_PROGRESS("Merging layouts via DLA", 100.0, NULL);

  /* Create the result */
  IGRAPH_CHECK(igraph_matrix_resize(res, allnodes, 2));
  respos=0;
  for (i=0; i<graphs; i++) {
    long int size=igraph_matrix_nrow(VECTOR(*coords)[i]);
    igraph_real_t xx=VECTOR(x)[i];
    igraph_real_t yy=VECTOR(y)[i];
    igraph_real_t rr=VECTOR(r)[i]/VECTOR(nr)[i];
    igraph_matrix_t *mat=VECTOR(*coords)[i];
    IGRAPH_ALLOW_INTERRUPTION();
    if (VECTOR(nr)[i]==0) { rr=1; }
    for (j=0; j<size; j++) {
      MATRIX(*res, respos, 0)=rr*(MATRIX(*mat, j, 0)-VECTOR(nx)[i]);
      MATRIX(*res, respos, 1)=rr*(MATRIX(*mat, j, 1)-VECTOR(ny)[i]);
      MATRIX(*res, respos, 0)+=xx;
      MATRIX(*res, respos, 1)+=yy;
      ++respos;
    }
  }

  RNG_END();
 
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

  /* The graph is not used, only its coordinates */
  IGRAPH_UNUSED(actg);

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

int igraph_i_layout_mds_step(igraph_real_t *to, const igraph_real_t *from,
			     int n, void *extra);

int igraph_i_layout_mds_single(const igraph_t* graph, igraph_matrix_t *res,
                               igraph_matrix_t *dist, long int dim);

int igraph_i_layout_mds_step(igraph_real_t *to, const igraph_real_t *from,
    int n, void *extra) {
  igraph_matrix_t* matrix = (igraph_matrix_t*)extra;
  IGRAPH_UNUSED(n);
  igraph_blas_dgemv_array(0, 1, matrix, from, 0, to);
  return 0;
}

/* MDS layout for a connected graph, with no error checking on the
 * input parameters. The distance matrix will be modified in-place. */
int igraph_i_layout_mds_single(const igraph_t* graph, igraph_matrix_t *res,
                               igraph_matrix_t *dist, long int dim) {

  long int no_of_nodes=igraph_vcount(graph);
  long int nev = dim;
  igraph_matrix_t vectors;
  igraph_vector_t values, row_means;
  igraph_real_t grand_mean;
  long int i, j, k;
  igraph_eigen_which_t which;

  /* Handle the trivial cases */
  if (no_of_nodes == 1) {
    IGRAPH_CHECK(igraph_matrix_resize(res, 1, dim));
    igraph_matrix_fill(res, 0);
    return IGRAPH_SUCCESS;
  }
  if (no_of_nodes == 2) {
    IGRAPH_CHECK(igraph_matrix_resize(res, 2, dim));
    igraph_matrix_fill(res, 0);
    for (j = 0; j < dim; j++)
      MATRIX(*res, 1, j) = 1;
    return IGRAPH_SUCCESS;
  }

  /* Initialize some stuff */
  IGRAPH_VECTOR_INIT_FINALLY(&values, no_of_nodes);
  IGRAPH_CHECK(igraph_matrix_init(&vectors, no_of_nodes, dim));
  IGRAPH_FINALLY(igraph_matrix_destroy, &vectors);

  /* Take the square of the distance matrix */
  for (i = 0; i < no_of_nodes; i++) {
    for (j = 0; j < no_of_nodes; j++) {
      MATRIX(*dist, i, j) *= MATRIX(*dist, i, j);
    }
  }

  /* Double centering of the distance matrix */
  IGRAPH_VECTOR_INIT_FINALLY(&row_means, no_of_nodes);
  igraph_vector_fill(&values, 1.0 / no_of_nodes);
  igraph_blas_dgemv(0, 1, dist, &values, 0, &row_means);
  grand_mean = igraph_vector_sum(&row_means) / no_of_nodes;
  igraph_matrix_add_constant(dist, grand_mean);
  for (i = 0; i < no_of_nodes; i++) {
    for (j = 0; j < no_of_nodes; j++) {
      MATRIX(*dist, i, j) -= VECTOR(row_means)[i] + VECTOR(row_means)[j];
      MATRIX(*dist, i, j) *= -0.5;
    }
  }
  igraph_vector_destroy(&row_means);
  IGRAPH_FINALLY_CLEAN(1);

  /* Calculate the top `dim` eigenvectors. */
  which.pos = IGRAPH_EIGEN_LA;
  which.howmany = (int) nev;
  IGRAPH_CHECK(igraph_eigen_matrix_symmetric(/*A=*/ 0, /*sA=*/ 0, 
			       /*fun=*/ igraph_i_layout_mds_step,
			       /*n=*/ (int) no_of_nodes, /*extra=*/ dist, 
			       /*algorithm=*/ IGRAPH_EIGEN_LAPACK,
			       &which, /*options=*/ 0, /*storage=*/ 0,
			       &values, &vectors));

  /* Calculate and normalize the final coordinates */
  for (j = 0; j < nev; j++) {
    VECTOR(values)[j] = sqrt(fabs(VECTOR(values)[j]));
  }
  IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, dim));
  for (i = 0; i < no_of_nodes; i++) {
    for (j = 0, k=nev-1; j < nev; j++, k--) {
      MATRIX(*res, i, k) = VECTOR(values)[j] * MATRIX(vectors, i, j);
    }
  }

  igraph_matrix_destroy(&vectors);
  igraph_vector_destroy(&values);
  IGRAPH_FINALLY_CLEAN(2);

  return IGRAPH_SUCCESS;
}

/**
 * \function igraph_layout_mds
 * \brief Place the vertices on a plane using multidimensional scaling.
 * 
 * </para><para>
 * This layout requires a distance matrix, where the intersection of
 * row i and column j specifies the desired distance between vertex i
 * and vertex j. The algorithm will try to place the vertices in a
 * space having a given number of dimensions in a way that approximates
 * the distance relations prescribed in the distance matrix. igraph
 * uses the classical multidimensional scaling by Torgerson; for more
 * details, see Cox &amp; Cox: Multidimensional Scaling (1994), Chapman
 * and Hall, London.
 *
 * </para><para>
 * If the input graph is disconnected, igraph will decompose it
 * first into its subgraphs, lay out the subgraphs one by one
 * using the appropriate submatrices of the distance matrix, and
 * then merge the layouts using \ref igraph_layout_merge_dla.
 * Since \ref igraph_layout_merge_dla works for 2D layouts only,
 * you cannot run the MDS layout on disconnected graphs for
 * more than two dimensions.
 * 
 * </para><para>
 * Warning: if the graph is symmetric to the exchange of two vertices
 * (as is the case with leaves of a tree connecting to the same parent),
 * classical multidimensional scaling may assign the same coordinates to
 * these vertices.
 *
 * \param graph A graph object.
 * \param res Pointer to an initialized matrix object. This will
 *        contain the result and will be resized if needed.
 * \param dist The distance matrix. It must be symmetric and this
 *        function does not check whether the matrix is indeed
 *        symmetric. Results are unspecified if you pass a non-symmetric
 *        matrix here. You can set this parameter to null; in this
 *        case, the shortest path lengths between vertices will be
 *        used as distances.
 * \param dim The number of dimensions in the embedding space. For
 *        2D layouts, supply 2 here.
 * \param options This argument is currently ignored, it was used for 
 *        ARPACK, but LAPACK is used now for calculating the eigenvectors.
 * \return Error code.
 * 
 * Added in version 0.6.
 * 
 * </para><para>
 * Time complexity: usually around O(|V|^2 dim).
 */

int igraph_layout_mds(const igraph_t* graph, igraph_matrix_t *res,
                      const igraph_matrix_t *dist, long int dim,
                      igraph_arpack_options_t *options) {
  long int i, no_of_nodes=igraph_vcount(graph);
  igraph_matrix_t m;
  igraph_bool_t conn;

  RNG_BEGIN();

  /* Check the distance matrix */
  if (dist && (igraph_matrix_nrow(dist) != no_of_nodes ||
      igraph_matrix_ncol(dist) != no_of_nodes)) {
    IGRAPH_ERROR("invalid distance matrix size", IGRAPH_EINVAL);
  }

  /* Check the number of dimensions */
  if (dim <= 1) {
    IGRAPH_ERROR("dim must be positive", IGRAPH_EINVAL);
  }
  if (dim > no_of_nodes) {
    IGRAPH_ERROR("dim must be less than the number of nodes", IGRAPH_EINVAL);
  }

  /* Copy or obtain the distance matrix */
  if (dist == 0) {
    IGRAPH_CHECK(igraph_matrix_init(&m, no_of_nodes, no_of_nodes));
    IGRAPH_FINALLY(igraph_matrix_destroy, &m);
    IGRAPH_CHECK(igraph_shortest_paths(graph, &m,
          igraph_vss_all(), igraph_vss_all(), IGRAPH_ALL));
  } else {
    IGRAPH_CHECK(igraph_matrix_copy(&m, dist));
    IGRAPH_FINALLY(igraph_matrix_destroy, &m);
    /* Make sure that the diagonal contains zeroes only */
    for (i = 0; i < no_of_nodes; i++)
      MATRIX(m, i, i) = 0.0;
  }

  /* Check whether the graph is connected */
  IGRAPH_CHECK(igraph_is_connected(graph, &conn, IGRAPH_WEAK));
  if (conn) {
    /* Yes, it is, just do the MDS */
    IGRAPH_CHECK(igraph_i_layout_mds_single(graph, res, &m, dim));
  } else {
    /* The graph is not connected, lay out the components one by one */
    igraph_vector_ptr_t layouts;
    igraph_vector_t comp, vertex_order;
    igraph_t subgraph;
    igraph_matrix_t *layout;
    igraph_matrix_t dist_submatrix;
    igraph_bool_t *seen_vertices;
    long int j, n, processed_vertex_count = 0;

    IGRAPH_VECTOR_INIT_FINALLY(&comp, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&vertex_order, no_of_nodes);

    IGRAPH_CHECK(igraph_vector_ptr_init(&layouts, 0));
    IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &layouts);
    igraph_vector_ptr_set_item_destructor(&layouts, (igraph_finally_func_t*)igraph_matrix_destroy);

    IGRAPH_CHECK(igraph_matrix_init(&dist_submatrix, 0, 0));
    IGRAPH_FINALLY(igraph_matrix_destroy, &dist_submatrix);

    seen_vertices = igraph_Calloc(no_of_nodes, igraph_bool_t);
    if (seen_vertices == 0)
      IGRAPH_ERROR("cannot calculate MDS layout", IGRAPH_ENOMEM);
    IGRAPH_FINALLY(igraph_free, seen_vertices);

    for (i = 0; i < no_of_nodes; i++) {
      if (seen_vertices[i])
        continue;

      /* This is a vertex whose component we did not lay out so far */
      IGRAPH_CHECK(igraph_subcomponent(graph, &comp, i, IGRAPH_ALL));
      /* Take the subgraph */
      IGRAPH_CHECK(igraph_induced_subgraph(graph, &subgraph, igraph_vss_vector(&comp),
                                           IGRAPH_SUBGRAPH_AUTO));
      IGRAPH_FINALLY(igraph_destroy, &subgraph);
      /* Calculate the submatrix of the distances */
      IGRAPH_CHECK(igraph_matrix_select_rows_cols(&m, &dist_submatrix,
            &comp, &comp));
      /* Allocate a new matrix for storing the layout */
      layout = igraph_Calloc(1, igraph_matrix_t);
      if (layout == 0)
        IGRAPH_ERROR("cannot calculate MDS layout", IGRAPH_ENOMEM);
      IGRAPH_FINALLY(igraph_free, layout);
      IGRAPH_CHECK(igraph_matrix_init(layout, 0, 0));
      IGRAPH_FINALLY(igraph_matrix_destroy, layout);
      /* Lay out the subgraph */
      IGRAPH_CHECK(igraph_i_layout_mds_single(&subgraph, layout, &dist_submatrix, dim));
      /* Store the layout */
      IGRAPH_CHECK(igraph_vector_ptr_push_back(&layouts, layout));
      IGRAPH_FINALLY_CLEAN(2);  /* ownership of layout taken by layouts */
      /* Free the newly created subgraph */
      igraph_destroy(&subgraph);
      IGRAPH_FINALLY_CLEAN(1);
      /* Mark all the vertices in the component as visited */
      n = igraph_vector_size(&comp);
      for (j = 0; j < n; j++) {
        seen_vertices[(long int)VECTOR(comp)[j]] = 1;
        VECTOR(vertex_order)[(long int)VECTOR(comp)[j]] = processed_vertex_count++;
      }
    }
    /* Merge the layouts - reusing dist_submatrix here */
    IGRAPH_CHECK(igraph_layout_merge_dla(0, &layouts, &dist_submatrix));
    /* Reordering the rows of res to match the original graph */
    IGRAPH_CHECK(igraph_matrix_select_rows(&dist_submatrix, res, &vertex_order));

    igraph_free(seen_vertices);
    igraph_matrix_destroy(&dist_submatrix);
    igraph_vector_ptr_destroy_all(&layouts);
    igraph_vector_destroy(&vertex_order);
    igraph_vector_destroy(&comp);
    IGRAPH_FINALLY_CLEAN(5);
  }

  RNG_END();

  igraph_matrix_destroy(&m);
  IGRAPH_FINALLY_CLEAN(1);

  return IGRAPH_SUCCESS;
}

/**
 * \function igraph_layout_bipartite
 * Simple layout for bipartite graphs
 * 
 * The layout is created by first placing the vertices in two rows,
 * according to their types. Then the positions within the rows are
 * optimized to minimize edge crossings, by calling \ref
 * igraph_layout_sugiyama().
 * 
 * \param graph The input graph.
 * \param types A boolean vector containing ones and zeros, the vertex
 *     types. Its length must match the number of vertices in the graph.
 * \param res Pointer to an initialized matrix, the result, the x and
 *     y coordinates are stored here.
 * \param hgap The preferred minimum horizontal gap between vertices
 *     in the same layer (i.e. vertices of the same type).
 * \param vgap  The distance between layers.
 * \param maxiter Maximum number of iterations in the crossing
 *     minimization stage. 100 is a reasonable default; if you feel
 *     that you have too many edge crossings, increase this.
 * \return Error code.
 * 
 * \sa \ref igraph_layout_sugiyama().
 */

int igraph_layout_bipartite(const igraph_t *graph, 
			    const igraph_vector_bool_t *types,
			    igraph_matrix_t *res, igraph_real_t hgap, 
			    igraph_real_t vgap, long int maxiter) {

  long int i, no_of_nodes=igraph_vcount(graph);
  igraph_vector_t layers;

  if (igraph_vector_bool_size(types) != no_of_nodes) { 
    IGRAPH_ERROR("Invalid vertex type vector size", IGRAPH_EINVAL);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&layers, no_of_nodes);
  for (i=0; i<no_of_nodes; i++) {
    VECTOR(layers)[i] = 1 - VECTOR(*types)[i];
  }

  IGRAPH_CHECK(igraph_layout_sugiyama(graph, res, /*extd_graph=*/ 0, 
				      /*extd_to_orig_eids=*/ 0, &layers, hgap,
				      vgap, maxiter, /*weights=*/ 0));
  
  igraph_vector_destroy(&layers);
  IGRAPH_FINALLY_CLEAN(1);
  
  return 0;
}
