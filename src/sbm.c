/* -*- mode: C -*-  */
/* vim:set ts=8 sw=2 sts=2 et: */
/* 
   IGraph R library.
   Copyright (C) 2003-2013  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_interface.h"
#include "igraph_vector.h"
#include "igraph_matrix.h"
#include "igraph_random.h"
#include "igraph_constructors.h"
#include "igraph_games.h"

#include <float.h>		/* for DBL_EPSILON */
#include <math.h> 		/* for sqrt */

/**
 * \function igraph_sbm_game
 * Sample from a stochastic block model
 *
 * This function samples graphs from a stochastic block
 * model by (doing the equivalent of) Bernoulli
 * trials for each potential edge with the probabilities
 * given by the Bernoulli rate matrix, \p pref_matrix.
 * See Faust, K., &amp; Wasserman, S. (1992a). Blockmodels:
 * Interpretation and evaluation. Social Networks, 14, 5-–61.
 *
 * </para><para>
 * The order of the vertex ids in the generated graph corresponds to
 * the \p block_sizes argument.
 *
 * \param graph The output graph.
 * \param n Number of vertices.
 * \param pref_matrix The matrix giving the Bernoulli rates.
 *     This is a KxK matrix, where K is the number of groups.
 *     The probability of creating an edge between vertices from
 *     groups i and j is given by element (i,j).
 * \param block_sizes An integer vector giving the number of
 *     vertices in each group.
 * \param directed Boolean, whether to create a directed graph. If
 *     this argument is false, then \p pref_matrix must be symmetric.
 * \param loops Boolean, whether to create self-loops.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|+K^2), where |V| is the number of
 * vertices, |E| is the number of edges, and K is the number of
 * groups.
 *
 * \sa \ref igraph_erdos_renyi_game() for a simple Bernoulli graph.
 *
 */

int igraph_sbm_game(igraph_t *graph, igraph_integer_t n, 
		    const igraph_matrix_t *pref_matrix,
		    const igraph_vector_int_t *block_sizes,
		    igraph_bool_t directed, igraph_bool_t loops) {

  int no_blocks=igraph_matrix_nrow(pref_matrix);
  int from, to, fromoff=0;
  igraph_real_t minp, maxp;
  igraph_vector_t edges;
  
  /* ------------------------------------------------------------ */
  /* Check arguments                                              */
  /* ------------------------------------------------------------ */

  if (igraph_matrix_ncol(pref_matrix) != no_blocks) {
    IGRAPH_ERROR("Preference matrix is not square", 
		 IGRAPH_NONSQUARE);
  }

  igraph_matrix_minmax(pref_matrix, &minp, &maxp);
  if (minp < 0 || maxp > 1) { 
    IGRAPH_ERROR("Connection probabilities must in [0,1]", IGRAPH_EINVAL);
  }

  if (n < 0) { 
    IGRAPH_ERROR("Number of vertices must be non-negative", IGRAPH_EINVAL);
  }

  if (!directed && !igraph_matrix_is_symmetric(pref_matrix)) {
    IGRAPH_ERROR("Preference matrix must be symmetric for undirected graphs",
		 IGRAPH_EINVAL);
  }

  if (igraph_vector_int_size(block_sizes) != no_blocks) {
    IGRAPH_ERROR("Invalid block size vector length", IGRAPH_EINVAL);
  }

  if (igraph_vector_int_min(block_sizes) < 0) {
    IGRAPH_ERROR("Block size must be non-negative", IGRAPH_EINVAL);
  }

  if (igraph_vector_int_sum(block_sizes) != n) {
    IGRAPH_ERROR("Block sizes must sum up to number of vertices", 
		 IGRAPH_EINVAL);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);  
  
  RNG_BEGIN();

  for (from = 0; from < no_blocks; from++) {
    int fromsize = VECTOR(*block_sizes)[from];
    int start = directed ? 0 : from;
    int i, tooff=0;
    for (i=0; i<start; i++) {
      tooff += VECTOR(*block_sizes)[i];
    }
    for (to = start; to < no_blocks; to++) {
      int tosize = VECTOR(*block_sizes)[to];
      igraph_real_t prob=MATRIX(*pref_matrix, from, to);
      double maxedges, last=RNG_GEOM(prob);
      if (directed && loops) {
	maxedges = fromsize * tosize;
	while (last < maxedges) {
	  int vto=floor(last/fromsize);
	  int vfrom=last - (igraph_real_t)vto * fromsize;
	  igraph_vector_push_back(&edges, fromoff + vfrom);
	  igraph_vector_push_back(&edges, tooff + vto);
	  last += RNG_GEOM(prob);
	  last += 1;
	}
      } else if (directed && !loops && from!=to) {
	maxedges = fromsize * tosize;
	while (last < maxedges) {
	  int vto=floor(last/fromsize);
	  int vfrom=last - (igraph_real_t)vto * fromsize;
	  igraph_vector_push_back(&edges, fromoff + vfrom);
	  igraph_vector_push_back(&edges, tooff + vto);
	  last += RNG_GEOM(prob);
	  last += 1;
	}	
      } else if (directed && !loops && from==to) {
	maxedges = fromsize * (fromsize-1);
	while (last < maxedges) {
	  int vto=floor(last/fromsize);
	  int vfrom=last - (igraph_real_t)vto * fromsize;
	  if (vfrom == vto) { vto=fromsize-1; }
	  igraph_vector_push_back(&edges, fromoff + vfrom);
	  igraph_vector_push_back(&edges, tooff + vto);
	  last += RNG_GEOM(prob);
	  last += 1;
	}
      } else if (!directed && loops && from!=to) {
	maxedges = fromsize * tosize;	
	while (last < maxedges) {
	  int vto=floor(last/fromsize);
	  int vfrom=last - (igraph_real_t)vto * fromsize;
	  igraph_vector_push_back(&edges, fromoff + vfrom);
	  igraph_vector_push_back(&edges, tooff + vto);
	  last += RNG_GEOM(prob);
	  last += 1;
	}
      } else if (!directed && loops && from==to) {
	maxedges = fromsize * (fromsize+1) / 2.0;
	while (last < maxedges) {
	  long int vto=floor((sqrt(8*last+1)-1)/2);
	  long int vfrom=last-(((igraph_real_t)vto)*(vto+1))/2;
	  igraph_vector_push_back(&edges, fromoff + vfrom);
	  igraph_vector_push_back(&edges, tooff + vto);
	  last += RNG_GEOM(prob);
	  last += 1;
	}	
      } else if (!directed && !loops && from!=to) {
	maxedges = fromsize * tosize;	
	while (last < maxedges) {
	  int vto=floor(last/fromsize);
	  int vfrom=last - (igraph_real_t)vto * fromsize;
	  igraph_vector_push_back(&edges, fromoff + vfrom);
	  igraph_vector_push_back(&edges, tooff + vto);
	  last += RNG_GEOM(prob);
	  last += 1;
	}
      } else /*!directed && !loops && from==to */ {
	maxedges = fromsize * (fromsize-1) / 2.0;
	while (last < maxedges) {
	  int vto=floor((sqrt(8*last+1)+1)/2);
	  int vfrom=last-(((igraph_real_t)vto)*(vto-1))/2;
	  igraph_vector_push_back(&edges, fromoff + vfrom);
	  igraph_vector_push_back(&edges, tooff + vto);
	  last += RNG_GEOM(prob);
	  last += 1;
	}
      }
      
      tooff += tosize;
    }
    fromoff += fromsize;
  }

  RNG_END();

  igraph_create(graph, &edges, n, directed);
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

/**
 * \function igraph_hsbm_game
 *
 */

int igraph_hsbm_game(igraph_t *graph, igraph_integer_t n, 
		     igraph_integer_t m, const igraph_vector_t *rho,
		     const igraph_matrix_t *C, igraph_real_t p) {
  
  int b, i, k=igraph_vector_size(rho);
  igraph_vector_t csizes;
  igraph_real_t sq_dbl_epsilon=sqrt(DBL_EPSILON);
  int no_blocks=n / m;
  igraph_vector_t edges;
  int offset=0;

  if (n < 1) { 
    IGRAPH_ERROR("`n' must be positive for HSBM", IGRAPH_EINVAL); 
  }
  if (m < 1) {
    IGRAPH_ERROR("`m' must be positive for HSBM", IGRAPH_EINVAL);
  }
  if ((long) n  % (long) m) {
    IGRAPH_ERROR("`n' must be a multiple of `m' for HSBM", IGRAPH_EINVAL);
  }
  if (!igraph_vector_isininterval(rho, 0, 1)) {
    IGRAPH_ERROR("`rho' must be between zero and one for HSBM", 
		 IGRAPH_EINVAL);
  }
  if (fabs(igraph_vector_sum(rho) - 1.0) > sq_dbl_epsilon) {
    IGRAPH_ERROR("`rho' must sum up to 1 for HSBM", IGRAPH_EINVAL);
  }
  if (igraph_matrix_nrow(C) != k || igraph_matrix_ncol(C) != k) {
    IGRAPH_ERROR("`C' dimensions must match `rho' dimensions in HSBM",
		 IGRAPH_EINVAL);
  }
  if (!igraph_matrix_is_symmetric(C)) {
    IGRAPH_ERROR("`C' must be a symmetric matrix", IGRAPH_EINVAL);
  }
  if (p < 0 || p > 1) {
    IGRAPH_ERROR("`p' must be a probability for HSBM", IGRAPH_EINVAL);
  }
  for (i=0; i<k; i++) {
    igraph_real_t s=VECTOR(*rho)[i] * m;
    if (fabs(floor(s)-s) > sq_dbl_epsilon) {
      IGRAPH_ERROR("`rho' * `m' is not integer in HSBM", IGRAPH_EINVAL);
    }
  }
  
  IGRAPH_VECTOR_INIT_FINALLY(&csizes, k);
  for (i=0; i<k; i++) {
    VECTOR(csizes)[i] = round(VECTOR(*rho)[i] * m);
  }
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

  RNG_BEGIN();

  /* Block models first */

  for (b=0; b<no_blocks; b++) {
    int from, to, fromoff=0;
    
    for (from = 0; from < k; from++) {
      int fromsize = VECTOR(csizes)[from];
      int i, tooff=0;
      for (i=0; i<from; i++) {
	tooff += VECTOR(csizes)[i];
      }
      for (to = from; to < k; to++) {
	int tosize = VECTOR(csizes)[to];
	igraph_real_t prob=MATRIX(*C, from, to);
	igraph_real_t maxedges;
	igraph_real_t last=RNG_GEOM(prob);
	if (from != to) {
	  maxedges = fromsize * tosize;	
	  while (last < maxedges) {
	    int vto=floor(last/fromsize);
	    int vfrom=last - (igraph_real_t)vto * fromsize;
	    igraph_vector_push_back(&edges, offset + fromoff + vfrom);
	    igraph_vector_push_back(&edges, offset + tooff + vto);
	    last += RNG_GEOM(prob);
	    last += 1;
	  }
	} else /* from==to */ {
	  maxedges = fromsize * (fromsize-1) / 2.0;
	  while (last < maxedges) {
	    int vto=floor((sqrt(8*last+1)+1)/2);
	    int vfrom=last-(((igraph_real_t)vto)*(vto-1))/2;
	    igraph_vector_push_back(&edges, offset + fromoff + vfrom);
	    igraph_vector_push_back(&edges, offset + tooff + vto);
	    last += RNG_GEOM(prob);
	    last += 1;
	  }
	}

	tooff += tosize;
      }
      fromoff += fromsize;
    }
    
    offset += m;
  }

  /* And now the rest, if not a special case */

  if (p == 1) {
    int fromoff=0, tooff=m;
    for (b=0; b<no_blocks; b++) {
      igraph_real_t fromsize = m;
      igraph_real_t tosize = n - tooff;
      int from, to;
      for (from=0; from<fromsize; from++) {
	for (to=0; to<tosize; to++) {
	  igraph_vector_push_back(&edges, fromoff + from);
	  igraph_vector_push_back(&edges, tooff + to);
	}
      }
      fromoff += m;
      tooff += m;
    }
  } else if (p > 0) {
    int fromoff=0, tooff=m;
    for (b=0; b<no_blocks; b++) {
      igraph_real_t fromsize = m;
      igraph_real_t tosize = n - tooff;
      igraph_real_t maxedges= fromsize * tosize;
      igraph_real_t last=RNG_GEOM(p);
      while (last < maxedges) {
	int vto = floor(last/fromsize);
	int vfrom = last - (igraph_real_t) vto * fromsize;
	igraph_vector_push_back(&edges, fromoff + vfrom);
	igraph_vector_push_back(&edges, tooff + vto);
	last += RNG_GEOM(p);
	last += 1;
      }
      
      fromoff += m;
      tooff += m;
    }
  }

  RNG_END();

  igraph_create(graph, &edges, n, /*directed=*/ 0);

  igraph_vector_destroy(&edges);
  igraph_vector_destroy(&csizes);
  IGRAPH_FINALLY_CLEAN(2);

  return 0;
}
