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

/* The original version of this file was written by Pascal Pons
   The original copyright notice follows here. The FSF address was
   fixed by Tamas Nepusz */

// File: walktrap.cpp
//-----------------------------------------------------------------------------
// Walktrap v0.2 -- Finds community structure of networks using random walks
// Copyright (C) 2004-2005 Pascal Pons
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
// 02110-1301 USA
//-----------------------------------------------------------------------------
// Author   : Pascal Pons 
// Email    : pons@liafa.jussieu.fr
// Web page : http://www.liafa.jussieu.fr/~pons/
// Location : Paris, France
// Time	    : June 2005
//-----------------------------------------------------------------------------
// see readme.txt for more details

#include "walktrap_graph.h"
#include "walktrap_communities.h"
#include <ctime>
#include <set>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include "igraph.h"

using namespace std;

/** 
 * \function igraph_community_walktrap
 * 
 * This function is the implementation of the Walktrap community
 * finding algorithm, see Pascal Pons, Matthieu Latapy: Computing
 * communities in large networks using random walks, 
 * http://arxiv.org/abs/physics/0512106
 * 
 * </para><para>
 * Currently the original C++ implementation is used in igraph, 
 * see http://www.liafa.jussieu.fr/~pons/index.php?item=prog&amp;item2=walktrap&amp;lang=en
 * I'm grateful to Matthieu Latapy and Pascal Pons for providing this 
 * source code.
 *
 * </para><para>
 * Note that the graph must not contain isolated vertices in order to
 * use this method.
 *
 * \param graph The input graph.
 * \param weights Numeric vector giving the weights of the edges. 
 *     If it is a NULL pointer then all edges will have equal
 *     weights. The weights are expected to be positive.
 * \param steps Integer constant, the length of the random walks.
 * \param merges Pointer to a matrix, the merges performed by the
 *     algorithm will be stored here (if not NULL). Each merge is a
 *     row in a two-column matrix and contains the ids of the merged
 *     clusters. Clusters are numbered from zero and cluster number 
 *     smaller than the number of nodes in the network belong to the
 *     individual vertices as singleton clusters. In each step a new
 *     cluster is created from two other clusters and its id will be 
 *     one larger than the largest cluster id so far. This means that 
 *     before the first merge we have \c n clusters (the number of
 *     vertices in the graph) numbered from zero to \c n-1. The first
 *     merge created cluster \c n, the second cluster \c n+1, etc.
 * \param modularity Pointer to a vector. If not NULL then the
 *     modularity score of the current clustering is stored here after
 *     each merge operation. 
 * \return Error code.
 * 
 * \sa \ref igraph_community_spinglass(), \ref
 * igraph_community_edge_betweenness(). 
 * 
 * Time complexity: O(|E||V|^2) in the worst case, O(|V|^2 log|V|) typically, 
 * |V| is the number of vertices, |E| is the number of edges.
 */

int igraph_community_walktrap(const igraph_t *graph, 
			      const igraph_vector_t *weights,
			      int steps,
			      igraph_matrix_t *merges,
			      igraph_vector_t *modularity) {

  long int no_of_nodes=(long int)igraph_vcount(graph);
  int length=steps;
  long max_memory=-1;
  
  Graph* G = new Graph;
  if (G->convert_from_igraph(graph, weights))
      IGRAPH_ERROR("isolated vertex found in graph", IGRAPH_EINVAL);
  
  if (merges) {
    IGRAPH_CHECK(igraph_matrix_resize(merges, no_of_nodes-1, 2));
  }
  if (modularity) {
    IGRAPH_CHECK(igraph_vector_resize(modularity, no_of_nodes));
	igraph_vector_null(modularity);
  }
  Communities C(G, length, max_memory, merges, modularity);
  
  while (!C.H->is_empty()) {
    IGRAPH_ALLOW_INTERRUPTION();
    C.merge_nearest_communities();
  }
  
  delete G;
  
  return 0;
}
