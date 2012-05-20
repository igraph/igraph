/* -*- mode: C -*-  */
/* vim:set ts=2 sw=2 sts=2 et: */
/* 
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
#include "igraph_structural.h"
#include "igraph_community.h"
#include "igraph_error.h"
#include "igraph_glpk_support.h"
#include "igraph_interrupt_internal.h"
#include "config.h"

#ifdef HAVE_GLPK
#include <glpk.h>
#endif

/**
 * \function igraph_community_optimal_modularity
 * Calculate the community structure with the highest modularity value
 * 
 * This function calculates the optimal community structure for a
 * graph, in terms of maximal modularity score. 
 * 
 * </para><para>
 * The calculation is done by transforming the modularity maximization
 * into an integer programming problem, and then calling the GLPK
 * library to solve that. Please see Ulrik Brandes et al.: On
 * Modularity Clustering, IEEE Transactions on Knowledge and Data
 * Engineering 20(2):172-188, 2008.
 * 
 * </para><para>
 * Note that modularity optimization is an NP-complete problem, and
 * all known algorithms for it have exponential time complexity. This
 * means that you probably don't want to run this function on larger
 * graphs. Graphs with up to fifty vertices should be fine, graphs
 * with a couple of hundred vertices might be possible.
 * 
 * \param graph The input graph. It is always treated as undirected.
 * \param modularity Pointer to a real number, or a null pointer.
 *        If it is not a null pointer, then a optimal modularity value 
 *        is returned here.
 * \param membership Pointer to a vector, or a null pointer. If not a
 *        null pointer, then the membership vector of the optimal
 *        community structure is stored here.
 * \return Error code.
 * 
 * \sa \ref igraph_modularity(), \ref igraph_community_fastgreedy()
 * for an algorithm that finds a local optimum in a greedy way. 
 * 
 * Time complexity: exponential in the number of vertices.
 * 
 * \example examples/simple/igraph_community_optimal_modularity.c
 */

int igraph_community_optimal_modularity(const igraph_t *graph,
					igraph_real_t *modularity,
					igraph_vector_t *membership,
					igraph_bool_t verbose) {
#ifndef HAVE_GLPK
  IGRAPH_ERROR("GLPK is not available", 
	       IGRAPH_UNIMPLEMENTED);    
#else

  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_bool_t directed=igraph_is_directed(graph);
  long int no_of_variables=no_of_nodes * (no_of_nodes+1)/2;
  long int i, j, k, st;
  int idx[] = { 0, 0, 0, 0 };
  double coef[] = { 0.0, 1.0, 1.0, -2.0 };

  igraph_vector_t degree;

  glp_prob *ip;
  glp_iocp parm;

  /* Special case */
  if (no_of_edges == 0) {
    if (modularity) {
      *modularity=IGRAPH_NAN;
    }
    if (membership) {
      IGRAPH_CHECK(igraph_vector_resize(membership, no_of_nodes));
      igraph_vector_null(membership);
    }
  }
  
  IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
  IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(), 
			     IGRAPH_ALL, IGRAPH_LOOPS));

  glp_term_out(GLP_OFF);
  ip = glp_create_prob();
  IGRAPH_FINALLY(glp_delete_prob, ip);

  glp_set_obj_dir(ip, GLP_MAX);
  st=glp_add_cols(ip, no_of_variables);

  /* variables are binary */
  for (i=0; i<no_of_variables; i++) {
    glp_set_col_kind(ip, st+i, GLP_BV);
  }

#define IDX(a,b) ((b)*((b)+1)/2+(a))

  /* reflexivity */
  for (i=0; i<no_of_nodes; i++) {
    glp_set_col_bnds(ip, st+IDX(i,i), GLP_FX, 1.0, 1.0);
  }
  
  /* transitivity */
  for (i=0; i<no_of_nodes; i++) {
    for (j=i+1; j<no_of_nodes; j++) {

      IGRAPH_ALLOW_INTERRUPTION();

      for (k=j+1; k<no_of_nodes; k++) {
	long int newrow=glp_add_rows(ip, 3);

	glp_set_row_bnds(ip, newrow, GLP_UP, 0.0, 1.0);
	idx[1]  = st+IDX(i,j); idx[2] = st+IDX(j,k); idx[3] = st+IDX(i,k);
	glp_set_mat_row(ip, newrow, 3, idx, coef);

	glp_set_row_bnds(ip, newrow+1, GLP_UP, 0.0, 1.0);
	idx[1] = st+IDX(i,j); idx[2] = st+IDX(i,k); idx[3] = st+IDX(j,k);
	glp_set_mat_row(ip, newrow+1, 3, idx, coef);

	glp_set_row_bnds(ip, newrow+2, GLP_UP, 0.0, 1.0);
	idx[1] = st+IDX(i,k); idx[2] = st+IDX(j,k); idx[3] = st+IDX(i,j);
	glp_set_mat_row(ip, newrow+2, 3, idx, coef);

      }
    }
  }

  /* objective function */
  for (i=0; i<no_of_nodes; i++) {
    for (j=i; j<no_of_nodes; j++) {
      igraph_real_t c=1.0/2.0/no_of_edges;
      igraph_bool_t con;
      igraph_real_t e;
      igraph_are_connected(graph, i, j, &con);
      if (!con && directed) { igraph_are_connected(graph, j, i, &con); }
      e= con ? 1.0 : 0.0;
			if (i == j) e *= 2;
      c *= e-VECTOR(degree)[i]*VECTOR(degree)[j] / 2.0 / no_of_edges;
      if (i != j) { c *= 2.0; }
      glp_set_obj_coef(ip, st+IDX(i,j), c);
    }
  }

  /* solve it */
  glp_init_iocp(&parm);
  parm.br_tech = GLP_BR_DTH;
  parm.bt_tech = GLP_BT_BLB;
  parm.presolve = GLP_ON;
  parm.binarize = GLP_ON;
  parm.cb_func = igraph_i_glpk_interruption_hook;
  IGRAPH_GLPK_CHECK(glp_intopt(ip, &parm), "Modularity optimization failed");
  
  /* store the results */
  if (modularity) {
    *modularity=glp_mip_obj_val(ip);
  }
  
  if (membership) {
    long int comm=0;	 /* id of the last community that was found */
    IGRAPH_CHECK(igraph_vector_resize(membership, no_of_nodes));
    for (i=0; i<no_of_nodes; i++) {
      
      IGRAPH_ALLOW_INTERRUPTION();
      
      for (j=0; j<i; j++) {
	int val=glp_mip_col_val(ip, st+IDX(j,i));
	if (val==1) {
	  VECTOR(*membership)[i]=VECTOR(*membership)[j];
	  break;
	}
      }
      if (j==i) {		/* new community */
	VECTOR(*membership)[i]=comm++;
      }
    }
  }
  
#undef IDX

  igraph_vector_destroy(&degree);
  glp_delete_prob(ip);
  IGRAPH_FINALLY_CLEAN(2);

  return 0;
  
#endif

}

