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
  igraph_vector_t deg;
  
  igraph_rng_seed(igraph_rng_default(), 42);

  igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, 100000, 
			  10000000, IGRAPH_DIRECTED, IGRAPH_LOOPS);
  
  igraph_vector_init(&deg, 0);

  BENCH("Degree, out, loops:      ", 25,
	igraph_degree(&graph, &deg, igraph_vss_all(), IGRAPH_OUT, 
		      IGRAPH_LOOPS);
	);

  BENCH("Degree, total, loops:    ", 25,
	igraph_degree(&graph, &deg, igraph_vss_all(), IGRAPH_ALL, 
		      IGRAPH_LOOPS);
	);

  BENCH("Degree, out, no loops:   ", 25,
	igraph_degree(&graph, &deg, igraph_vss_all(), IGRAPH_OUT, 
		      IGRAPH_NO_LOOPS);
	);

  BENCH("Degree, total, no loops: ", 25,
	igraph_degree(&graph, &deg, igraph_vss_all(), IGRAPH_ALL, 
		      IGRAPH_NO_LOOPS);
	);
  
  igraph_vector_destroy(&deg);
  igraph_destroy(&graph);
    
  return 0;
}
