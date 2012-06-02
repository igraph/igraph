/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
  
  igraph_t g, g2;
  igraph_bool_t iso;
  
  igraph_lcf(&g, 12, 5, -5, 6, 0);
  igraph_famous(&g2, "franklin");
  
  igraph_isomorphic_vf2(&g, &g2, 
			/*vertex.color1=*/ 0, /*vertex.color2=*/ 0, 
			/*edge.color1=*/ 0, /*edge.color2=*/ 0,
			&iso, 0, 0, 0, 0, 0);
  if (!iso) {
    printf("OOOPS!\n");
    return 1;
  }
  
  igraph_destroy(&g);
  igraph_destroy(&g2);
  
  return 0;
}
