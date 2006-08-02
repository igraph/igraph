/* -*- mode: C -*-  */
/* 
   IGraph library.
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include <igraph.h>

int main() {
  
  igraph_vector_t v;
  igraph_t g;
  
  igraph_vector_init_int(&v, 2, 5, 5);
  igraph_lattice(&g, &v, 1, IGRAPH_DIRECTED, 1 /*mutual*/, 0 /*circular*/);
  igraph_to_undirected(&g, IGRAPH_TO_UNDIRECTED_COLLAPSE);
  igraph_write_graph_edgelist(&g, stdout);

  igraph_destroy(&g);
  igraph_vector_destroy(&v);
  
  return 0;
}
