/* -*- mode: C -*-  */
/* 
   IGraph library.
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

int main() {
  
  igraph_t g;
  igraph_vector_t v;
  long int i;
  igraph_integer_t value, retcode;
  
  igraph_star(&g, 100, IGRAPH_STAR_UNDIRECTED, 0);
  
  igraph_vector_init(&v, 0);
  igraph_eigenvector_centrality(&g, &v, &value, &retcode, /*vmult*/0, /*aupdate*/0, 
				/*norm=*/0, /*tol=*/0, /*maxiter=*/300, /*ncv=*/3);

  if (retcode != 0) {
    return 1;
  }

  for (i=0; i<igraph_vector_size(&v); i++) {
    printf(" %.3f", fabs(VECTOR(v)[i]));
  }
  printf("\n");
  
  igraph_vector_destroy(&v);
  igraph_destroy(&g);
  
  return 0;
}
