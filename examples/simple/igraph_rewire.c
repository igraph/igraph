/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

int main() {
  
  igraph_t g;
  igraph_vector_t res;
  long i, n;
  
  igraph_tree(&g, 10, 3, IGRAPH_TREE_OUT);
  igraph_rewire(&g, 1000, IGRAPH_REWIRING_SIMPLE);
  
  n=igraph_vcount(&g);
  igraph_vector_init(&res, 0);
  igraph_degree(&g, &res, igraph_vss_all(), IGRAPH_IN, 0);
  for (i=0; i<n; i++)
    printf("%ld ", (long)igraph_vector_e(&res, i));
  printf("\n");
  
  igraph_degree(&g, &res, igraph_vss_all(), IGRAPH_OUT, 0);
  for (i=0; i<n; i++)
    printf("%ld ", (long)igraph_vector_e(&res, i));
  printf("\n");
    
  igraph_destroy(&g);
  igraph_vector_destroy(&res);
  
  return 0;
}
