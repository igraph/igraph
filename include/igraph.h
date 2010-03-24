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

#ifndef IGRAPH_H
#define IGRAPH_H

#ifndef _GNU_SOURCE
# define _GNU_SOURCE 1
#endif

#include "igraph_version.h"
#include "igraph_memory.h"
#include "igraph_error.h"
#include "igraph_random.h"
#include "igraph_interrupt.h"
#include "igraph_progress.h"

#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_matrix.h"
#include "igraph_array.h"
#include "igraph_dqueue.h"
#include "igraph_stack.h"
#include "igraph_heap.h"
#include "igraph_psumtree.h"
#include "igraph_strvector.h"
#include "igraph_vector_ptr.h"
#include "igraph_spmatrix.h"

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
#include "igraph_attributes.h"
#include "igraph_arpack.h"
#include "igraph_mixing.h"

#endif
