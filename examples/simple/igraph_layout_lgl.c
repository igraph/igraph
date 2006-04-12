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
#include <math.h>

int main() {

  igraph_t g;
  igraph_matrix_t coords;
  igraph_real_t vc;
  
  igraph_tree(&g, 100, 3, IGRAPH_TREE_UNDIRECTED);
/*   igraph_barabasi_game(&g, 1000, 1, 0, 0, IGRAPH_UNDIRECTED); */
  igraph_matrix_init(&coords, 0, 0);
  vc=igraph_vcount(&g);
  igraph_layout_lgl(&g, &coords, 
		    /* maxiter */    150, 
		    /* maxdelta */   vc,
		    /* area */       vc*vc,
		    /* coolexp */    1.5,
		    /* repulserad */ vc*vc*vc,
		    /* cellsize */   sqrt(sqrt(vc)),
		    /* root */       0);
  
  igraph_matrix_destroy(&coords);
  igraph_destroy(&g);
  return 0;
}
