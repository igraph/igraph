/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2014  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA
   
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
#include <stdio.h>

int main() {

  IGRAPH(g,  IGRAPH_DIRECTED,   0,1, 1,2, 2,3, 3,0);
  IGRAPH(g2, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,3, 3,0);

  igraph_write_graph_edgelist(&g, stdout);
  printf("--\n");
  igraph_write_graph_edgelist(&g2, stdout);
  
  igraph_destroy(&g);
  igraph_destroy(&g2);
  
  return 0;
}
