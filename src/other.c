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
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/

#include "igraph.h"
#include "types.h"

/**
 * \ingroup nongraph
 * \function igraph_running_mean
 * \brief Calculates the running mean of a vector
 * 
 * The running mean is defined by the mean of the
 * previous <parameter>binwidth</parameter> values.
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

int igraph_running_mean(const vector_t *data, vector_t *res, 
			integer_t binwidth) {

  double sum=0;
  long int i;

  /* Check */
  if (vector_size(data) < binwidth) {
    IGRAPH_ERROR("Vector too short for this binwidth", IGRAPH_EINVAL); 
  }

  /* Memory for result */

  IGRAPH_CHECK(vector_resize(res, (long int)(vector_size(data)-binwidth+1)));
  
  /* Initial bin */
  for (i=0; i<binwidth; i++) {
    sum += VECTOR(*data)[i];
  }
  
  VECTOR(*res)[0]=sum/binwidth;
  
  for (i=1; i<vector_size(data)-binwidth+1; i++) {
    sum -= VECTOR(*data)[i-1];
    sum += VECTOR(*data)[ (long int)(i+binwidth-1)];
    VECTOR(*res)[i] = sum/binwidth;
  }
  
  return 0;
}


