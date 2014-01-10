/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2014  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include "igraph_games.h"
#include "igraph_random.h"
#include "igraph_constructors.h"
#include "igraph_lapack.h"

int igraph_dot_product_game(igraph_t *graph, const igraph_matrix_t *vecs,
			    igraph_bool_t directed) {

  igraph_integer_t nrow=igraph_matrix_nrow(vecs);
  igraph_integer_t ncol=igraph_matrix_ncol(vecs);
  int i, j;
  igraph_vector_t edges;
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    
  RNG_BEGIN();

  for (i = 0; i < ncol; i++) {
    int from=directed ? 0 : i+1;
    igraph_vector_t v1;
    igraph_vector_view(&v1, &MATRIX(*vecs, 0, i), nrow);
    for (j = from; j < ncol; j++) {
      if (i==j) { continue; }
      igraph_real_t prob;
      igraph_vector_t v2;
      igraph_vector_view(&v2, &MATRIX(*vecs, 0, j), nrow);
      igraph_lapack_ddot(&v1, &v2, &prob);
      if (prob < 0) { 
	IGRAPH_WARNING("Negative connection probability in "
		       "dot-product graph");
      } else if (prob > 1) {
	IGRAPH_WARNING("Greater than 1 connection probability in "
		       "dot-product graph");
	IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	IGRAPH_CHECK(igraph_vector_push_back(&edges, j));
      } else if (RNG_UNIF01() < prob) { 
	IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
	IGRAPH_CHECK(igraph_vector_push_back(&edges, j));
      }
    }
  }

  RNG_END();
  
  igraph_create(graph, &edges, ncol, directed);
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}
