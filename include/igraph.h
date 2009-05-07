/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2003, 2004, 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#ifndef RESTGAME_H
#define RESTGAME_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

#ifndef _GNU_SOURCE
# define _GNU_SOURCE
#endif

#include "igraph_types.h"
#include "igraph_error.h"
#include "igraph_interrupt.h"
#include "igraph_arpack.h"

#include <stdio.h> 		/* FILE */

#include "igraph_constants.h"
#include "igraph_datatype.h"
#include "igraph_iterators.h"
#include "igraph_interface.h"
#include "igraph_constructors.h"
#include "igraph_games.h"
#include "igraph_centrality.h"
#include "igraph_paths.h"
#include "igraph_components.h"
#include "igraph_structural.h"
#include "igraph_transitivity.h"
#include "igraph_neighborhood.h"
#include "igraph_topology.h"
#include "igraph_bipartite.h"
#include "igraph_cliques.h"
#include "igraph_layout.h"
#include "igraph_visitor.h"
#include "igraph_community.h"
#include "igraph_conversion.h"
#include "igraph_foreign.h"
#include "igraph_motifs.h"
#include "igraph_operators.h"
#include "igraph_flow.h"
#include "igraph_revolver.h"
#include "igraph_nongraph.h"
#include "igraph_cocitation.h"
#include "igraph_adjlist.h"
#include "igraph_progress.h"

/* -------------------------------------------------- */
/* Eigenvectors and eigenvalues                       */
/* -------------------------------------------------- */

int igraph_eigen_tred2(const igraph_matrix_t *A,
		       igraph_vector_t *D,
		       igraph_vector_t *E,
		       igraph_matrix_t *Z);

int igraph_eigen_tql2(igraph_vector_t *D,
		      igraph_vector_t *E,
		      igraph_matrix_t *Z);

int igraph_eigen_tred1(const igraph_matrix_t *A,
		       igraph_vector_t *D,
		       igraph_vector_t *E2);

int igraph_eigen_tqlrat(igraph_vector_t *D,
			igraph_vector_t *E2);

int igraph_eigen_rs(const igraph_matrix_t *A,
		    igraph_vector_t *values,
		    igraph_matrix_t *vectors);

#include "igraph_attributes.h"

__END_DECLS
  
#endif
