/* -*- mode: C -*-  */
/* 
   IGraph library.
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
#include "types.h"
#include "memory.h"
#include "config.h"
#include <math.h>
#include <stdarg.h>

/**
 * \ingroup nongraph
 * \function igraph_running_mean
 * \brief Calculates the running mean of a vector.
 * 
 * </para><para>
 * The running mean is defined by the mean of the
 * previous \p binwidth values.
 * \param data The vector containing the data.
 * \param res The vector containing the result. This should be
 *        initialized before calling this function and will be
 *        resized. 
 * \param binwidth Integer giving the width of the bin for the running
 *        mean calculation.
 * \return Error code.
 * 
 * Time complexity: O(n),
 * n is the length of
 * the data vector.
 */

int igraph_running_mean(const igraph_vector_t *data, igraph_vector_t *res, 
			igraph_integer_t binwidth) {

  double sum=0;
  long int i;

  /* Check */
  if (igraph_vector_size(data) < binwidth) {
    IGRAPH_ERROR("Vector too short for this binwidth", IGRAPH_EINVAL); 
  }

  /* Memory for result */

  IGRAPH_CHECK(igraph_vector_resize(res, (long int)(igraph_vector_size(data)-binwidth+1)));
  
  /* Initial bin */
  for (i=0; i<binwidth; i++) {
    sum += VECTOR(*data)[i];
  }
  
  VECTOR(*res)[0]=sum/binwidth;
  
  for (i=1; i<igraph_vector_size(data)-binwidth+1; i++) {
    IGRAPH_ALLOW_INTERRUPTION();
    sum -= VECTOR(*data)[i-1];
    sum += VECTOR(*data)[ (long int)(i+binwidth-1)];
    VECTOR(*res)[i] = sum/binwidth;
  }
  
  return 0;
}


/**
 * \ingroup nongraph
 * \function igraph_convex_hull
 * \brief Determines the convex hull of a given set of points in the 2D plane
 *
 * </para><para>
 * The convex hull is determined by the Graham scan algorithm.
 * See the following reference for details:
 * 
 * </para><para>
 * Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and Clifford
 * Stein. Introduction to Algorithms, Second Edition. MIT Press and
 * McGraw-Hill, 2001. ISBN 0262032937. Pages 949-955 of section 33.3:
 * Finding the convex hull.
 * 
 * \param data vector containing the coordinates. The length of the
 *        vector must be even, since it contains X-Y coordinate pairs.
 * \param resverts the vector containing the result, e.g. the vector of
 *        vertex indices used as the corners of the convex hull. Supply
 *        \c NULL here if you are only interested in the coordinates of
 *        the convex hull corners.
 * \param rescoords the matrix containing the coordinates of the selected
 *        corner vertices. Supply \c NULL here if you are only interested in
 *        the vertex indices.
 * \return Error code:
 *         \c IGRAPH_ENOMEM: not enough memory
 * 
 * Time complexity: O(n log(n)) where n is the number of vertices
 */
int igraph_convex_hull(const igraph_matrix_t *data, igraph_vector_t *resverts,
		       igraph_matrix_t *rescoords) {
  igraph_integer_t no_of_nodes;
  long int i, pivot_idx=0, last_idx, before_last_idx, next_idx, j;
  igraph_real_t* angles;
  igraph_vector_t stack;
  igraph_indheap_t order;
  igraph_real_t px, py, cp;
  
  no_of_nodes=igraph_matrix_nrow(data);
  if (igraph_matrix_ncol(data) != 2) {
    IGRAPH_ERROR("matrix must have 2 columns", IGRAPH_EINVAL);
  }
  if (no_of_nodes == 0) {
    if (resverts != 0) {
      IGRAPH_CHECK(igraph_vector_resize(resverts, 0));
    } 
    if (rescoords != 0) {
      IGRAPH_CHECK(igraph_matrix_resize(rescoords, 0, 2));
    }
    /**************************** this is an exit here *********/
    return 0;
  }
    
  angles=igraph_Calloc(no_of_nodes, igraph_real_t);
  if (!angles) IGRAPH_ERROR("not enough memory for angle array", IGRAPH_ENOMEM);
  IGRAPH_FINALLY(free, angles);
  
  IGRAPH_VECTOR_INIT_FINALLY(&stack, 0);
  
  /* Search for the pivot vertex */
  for (i=1; i<no_of_nodes; i++) {
    if (MATRIX(*data, i, 1)<MATRIX(*data, pivot_idx, 1))
      pivot_idx=i;
    else if (MATRIX(*data, i, 1) == MATRIX(*data, pivot_idx, 1) &&
	     MATRIX(*data, i, 0) < MATRIX(*data, pivot_idx, 0))
      pivot_idx=i;
  }
  px=MATRIX(*data, pivot_idx, 0);
  py=MATRIX(*data, pivot_idx, 1);
  
  /* Create angle array */
  for (i=0; i<no_of_nodes; i++) {
    if (i == pivot_idx) {
      /* We can't calculate the angle of the pivot point with itself,
       * so we use 10 here. This way, after sorting the angle vector,
       * the pivot point will always be the first one, since the range
       * of atan2 is -3.14..3.14 */
      angles[i] = 10;
    } else {
      angles[i] = atan2(MATRIX(*data, i, 1)-py,
			MATRIX(*data, i, 0)-px);
    }
  }

  IGRAPH_CHECK(igraph_indheap_init_array(&order, angles, no_of_nodes));
  IGRAPH_FINALLY(igraph_indheap_destroy, &order);
  
  igraph_Free(angles);
  IGRAPH_FINALLY_CLEAN(1);

  if (no_of_nodes == 1) {
    IGRAPH_CHECK(igraph_vector_push_back(&stack, 0));
    igraph_indheap_delete_max(&order);
  } else {
    /* Do the trick */
    IGRAPH_CHECK(igraph_vector_push_back(&stack, igraph_indheap_max_index(&order)-1));
    igraph_indheap_delete_max(&order);
    IGRAPH_CHECK(igraph_vector_push_back(&stack, igraph_indheap_max_index(&order)-1));
    igraph_indheap_delete_max(&order);
    
    j=2;
    while (!igraph_indheap_empty(&order)) {
      /* Determine whether we are at a left or right turn */
      last_idx=VECTOR(stack)[j-1];
      before_last_idx=VECTOR(stack)[j-2];
      next_idx=(long)igraph_indheap_max_index(&order)-1;
      igraph_indheap_delete_max(&order);
      cp=(MATRIX(*data, last_idx, 0)-MATRIX(*data, before_last_idx, 0))*
	(MATRIX(*data, next_idx, 1)-MATRIX(*data, before_last_idx, 1))-
	(MATRIX(*data, next_idx, 0)-MATRIX(*data, before_last_idx, 0))*
	(MATRIX(*data, last_idx, 1)-MATRIX(*data, before_last_idx, 1));
      /*
       printf("B L N cp: %d, %d, %d, %f [", before_last_idx, last_idx, next_idx, (float)cp);
       for (k=0; k<j; k++) printf("%ld ", (long)VECTOR(stack)[k]);
       printf("]\n");
       */
      if (cp == 0) {
	/* The last three points are collinear. Replace the last one in
	 * the stack to the newest one */
	VECTOR(stack)[j-1]=next_idx;
      } else if (cp < 0) {
	/* We are turning into the right direction */
	IGRAPH_CHECK(igraph_vector_push_back(&stack, next_idx));
	j++;
      } else {
	/* No, skip back until we're okay */
	while (cp >= 0 && j > 2) {
	  igraph_vector_pop_back(&stack);
	  j--;
	  last_idx=VECTOR(stack)[j-1];
	  before_last_idx=VECTOR(stack)[j-2];
	  cp=(MATRIX(*data, last_idx, 0)-MATRIX(*data, before_last_idx, 0))*
	    (MATRIX(*data, next_idx, 1)-MATRIX(*data, before_last_idx, 1))-
	    (MATRIX(*data, next_idx, 0)-MATRIX(*data, before_last_idx, 0))*
	    (MATRIX(*data, last_idx, 1)-MATRIX(*data, before_last_idx, 1));
	}
	IGRAPH_CHECK(igraph_vector_push_back(&stack, next_idx));
	j++;
      }
    }
  }
  
  /* Create result vector */
  if (resverts != 0) {
    igraph_vector_clear(resverts);
    IGRAPH_CHECK(igraph_vector_append(resverts, &stack));
  } 
  if (rescoords != 0) {
    igraph_matrix_select_rows(data, rescoords, &stack);
  }
  
  /* Free everything */
  igraph_vector_destroy(&stack);
  igraph_indheap_destroy(&order);
  IGRAPH_FINALLY_CLEAN(2);
  
  return 0;
}


/**
 * Internal function, floating point division
 * Used only in compilers not supporting INFINITY and HUGE_VAL to create
 * infinity values
 */
double igraph_i_fdiv(const double a, const double b) 
{
   return a / b;
}

