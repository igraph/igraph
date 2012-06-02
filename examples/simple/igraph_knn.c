/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
  igraph_vector_t v, v2;
  
  igraph_vector_init(&v, 0);
  igraph_vector_init(&v2, 0);
  igraph_ring(&g, 10, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/ 1);
  igraph_avg_nearest_neighbor_degree(&g, igraph_vss_all(), &v, &v2, 
				     /*weights=*/ 0);
  
  igraph_destroy(&g);
  igraph_vector_destroy(&v);
  igraph_vector_destroy(&v2);
  
  return 0;
}

