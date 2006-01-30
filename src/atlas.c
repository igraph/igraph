/* -*- mode: C -*-  */
/* 
   IGraph R package.
   Copyright (C) 2006  Gabor Csardi <csardi@rmki.kfki.hu>
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
#include "atlas-edges.h"

int igraph_atlas(igraph_t *graph, int number) {
  
  long int pos, n, e;
  const igraph_vector_t v=IGRAPH_VECTOR_NULL;

  if (number < 0 ||
      number >= sizeof(igraph_i_atlas_edges_pos)/sizeof(long int)) {
    IGRAPH_ERROR("No such graph in atlas", IGRAPH_EINVAL);
  }

  pos=igraph_i_atlas_edges_pos[number];
  n=igraph_i_atlas_edges[pos];
  e=igraph_i_atlas_edges[pos+1];
  
  IGRAPH_CHECK(igraph_create(graph, 
			     igraph_vector_view(&v,igraph_i_atlas_edges+pos+2, 
						e*2),  
			     n, IGRAPH_UNDIRECTED));
  
  return 0;
}
