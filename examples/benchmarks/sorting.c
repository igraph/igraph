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

/* Benchmark radix sorting vs. quick sorting, on two test cases,
   one with randomly shuffled and one with ordered vectors. The
   latter (almost) happens when a couple of new edges are added
   to a graph. */

#define N 1000000

int main() {

  igraph_vector_t from, to, res;
  igraph_vector_init(&res, 0);

  igraph_rng_seed(igraph_rng_default(), 42);

  igraph_vector_init_seq(&from, 0, N-1);
  igraph_vector_init_seq(&to,   0, N-1);
  igraph_vector_shuffle(&from);
  igraph_vector_shuffle(&to);
  
  BENCH("1 Sorting, quick sort, random vector",
	igraph_vector_sort(&from); igraph_vector_sort(&to));

  igraph_vector_destroy(&from);
  igraph_vector_destroy(&to);

  /* ------------------------------------------------------------------ */
  
  igraph_rng_seed(igraph_rng_default(), 42);
  
  igraph_vector_init_seq(&from, 0, N-1);
  igraph_vector_init_seq(&to,   0, N-1);
  igraph_vector_shuffle(&from);
  igraph_vector_shuffle(&to);
  
  BENCH("2 Sorting, radix sort, random vector",
	igraph_vector_order(&from, &to, &res, N-1));
  
  igraph_vector_destroy(&from);
  igraph_vector_destroy(&to);
  
  /* ------------------------------------------------------------------ */

  igraph_vector_init_seq(&from, 0, N-1);
  igraph_vector_init_seq(&to,   0, N-1);
  
  BENCH("3 Sorting, quick sort, sorted vector",
	igraph_vector_sort(&from); igraph_vector_sort(&to));

  igraph_vector_destroy(&from);
  igraph_vector_destroy(&to);

  /* ------------------------------------------------------------------ */

  igraph_vector_init_seq(&from, 0, N-1);
  igraph_vector_init_seq(&to,   0, N-1);
  
  BENCH("4 Sorting, radix sort, sorted vector",
	igraph_vector_order(&from, &to, &res, N-1));
  
  igraph_vector_destroy(&from);
  igraph_vector_destroy(&to);

  /* ------------------------------------------------------------------ */
  
  igraph_vector_destroy(&res);
  
  return 0;
}
