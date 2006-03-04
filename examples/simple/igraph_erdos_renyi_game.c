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

  igraph_t g;
  int i, ret;
  
  /* G(n,p) */
  igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 10, 0.0, 
			  IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
  if (igraph_ecount(&g) != 0) {
    return 1;
  }
  if (igraph_is_directed(&g)) {
    return 2;
  }
  igraph_destroy(&g);
  
  igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 10, 1.0,
			  IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
  if (igraph_ecount(&g) != 10*9) {
    return 3;
  }
  if (!igraph_is_directed(&g)) {
    return 4;
  }
  igraph_destroy(&g);

  igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 10, 0.5,
			  IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
  igraph_destroy(&g);

  /* G(n,m) */
  

  return 0;
}
