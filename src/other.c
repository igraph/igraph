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

#include "types.h"

int igraph_running_mean(vector_t *data, vector_t *res, integer_t binwidth) {

  double sum=0;
  long int i;

  /* Check */
  if (vector_size(data) < binwidth) {
    igraph_error("Vector too short for this binwidth");
  }

  /* Memory for result */

  vector_resize(res, (long int)(vector_size(data)-binwidth+1));
  
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


