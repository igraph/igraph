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
   The original copyright notice follows here */

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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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

int igraph_walktrap_community(const igraph_t *graph, 
			      const igraph_vector_t *weights,
			      int steps,
			      igraph_matrix_t *merges,
			      igraph_vector_t *modularity) {

  long int no_of_nodes=(long int)igraph_vcount(graph);
  int length=steps;
  long max_memory=-1;
  
  Graph* G = new Graph;
  G->convert_from_igraph(graph, weights);
  
  if (merges) {
    IGRAPH_CHECK(igraph_matrix_resize(merges, no_of_nodes-1, 2));
  }
  if (modularity) {
    IGRAPH_CHECK(igraph_vector_resize(modularity, no_of_nodes));
  }
  Communities C(G, length, max_memory, merges, modularity);
  
  while (!C.H->is_empty()) {
    C.merge_nearest_communities();
  }
  
  delete G;
  
  return 0;
}
