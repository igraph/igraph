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
  igraph_real_t res;

  /* Trivial cases */

  igraph_ring(&g, 100, IGRAPH_UNDIRECTED, 0, 0);
  igraph_reciprocity(&g, &res, 0);
  igraph_destroy(&g);

  if (res != 1) {
    return 1;
  }

  /* Small test graph */

  igraph_small(&g, 0, IGRAPH_DIRECTED, 
	       0,  1,  0,  2,  0,  3,  1,  0,  2,  3,  3,  2, -1);
  
  igraph_reciprocity(&g, &res, 0);
  igraph_destroy(&g);

  if (res != 0.5) {
    fprintf(stderr, "%f != %f\n", res, 0.5);
    return 2;
  }

  return 0;
}
