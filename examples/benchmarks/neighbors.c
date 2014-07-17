/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2014  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA
   
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

#include <igraph.h>
#include "bench.h"

int main() {

  igraph_t graph;
  igraph_adjlist_t al;
  
  igraph_rng_seed(igraph_rng_default(), 42);

  igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 100000, 
			  10000000, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
  
  BENCH("Create an adjacency list, degree 100: ", 10,
	igraph_adjlist_init(&graph, &al, IGRAPH_ALL);
	igraph_adjlist_destroy(&al);
	);
  
  igraph_destroy(&graph);

  igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 100000, 
			  1000000, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

  BENCH("Create an adjacency list, degree 10 : ", 10,
	igraph_adjlist_init(&graph, &al, IGRAPH_ALL);
	igraph_adjlist_destroy(&al);
	);
  
  igraph_destroy(&graph);
    
  return 0;
}
