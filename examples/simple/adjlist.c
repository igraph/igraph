/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2008  Gabor Csardi <csardi@rmki.kfki.hu>
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

#include <igraph.h>

int main() {

  igraph_t g, g2;
  igraph_adjlist_t adjlist;
  igraph_bool_t iso;
  
  /* Directed */
  igraph_tree(&g, 42, 3, IGRAPH_TREE_OUT);
  igraph_adjlist_init(&g, &adjlist, IGRAPH_OUT);
  igraph_adjlist(&g2, &adjlist, /*directed=*/ 1, /*duplicate=*/ 0);
  igraph_isomorphic_bliss(&g, &g2, &iso, 0, 0, IGRAPH_BLISS_F, IGRAPH_BLISS_F, 0, 0);
  if (!iso) { return 1; }
  igraph_adjlist_destroy(&adjlist);
  igraph_destroy(&g2);
  igraph_destroy(&g);

  /* Undirected */
  igraph_tree(&g, 42, 3, IGRAPH_TREE_UNDIRECTED);
  igraph_adjlist_init(&g, &adjlist, IGRAPH_OUT);
  igraph_adjlist(&g2, &adjlist, /*directed=*/ 0, /*duplicate=*/ 1);
  igraph_isomorphic_bliss(&g, &g2, &iso, 0, 0, IGRAPH_BLISS_F, IGRAPH_BLISS_F, 0, 0);
  if (!iso) { return 1; }
  igraph_adjlist_destroy(&adjlist);
  igraph_destroy(&g2);
  igraph_destroy(&g);
			  
  return 0;
}

