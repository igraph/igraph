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
  long int i;
  

  /* DIRECTED */
  
  igraph_star(&g, 10, IGRAPH_STAR_OUT, 0);  

  for (i=0; i<100; i++) {
    igraph_es_t es;
    igraph_eit_t it;
    igraph_es_pairs_small(&es, IGRAPH_DIRECTED, 
			  0,1,0,2,0,5,0,2,0,3,0,4,0,7,0,9, -1);
    igraph_eit_create(&g, es, &it);
    IGRAPH_EIT_RESET(it);
    while (!IGRAPH_EIT_END(it)) {
      IGRAPH_EIT_GET(it);
      IGRAPH_EIT_NEXT(it);
    }
    igraph_eit_destroy(&it);
    igraph_es_destroy(&es);
  }

  igraph_destroy(&g);

  /* UNDIRECTED */
  
  igraph_star(&g, 10, IGRAPH_STAR_UNDIRECTED, 0);
  
  for (i=0; i<100; i++) {
    igraph_es_t es;
    igraph_eit_t it;
    igraph_es_pairs_small(&es, IGRAPH_DIRECTED,
			  0,1,2,0,5,0,0,2,3,0,0,4,7,0,0,9, -1);
    igraph_eit_create(&g, es, &it);
    IGRAPH_EIT_RESET(it);
    while (!IGRAPH_EIT_END(it)) {
      IGRAPH_EIT_GET(it);
      IGRAPH_EIT_NEXT(it);
    }
    igraph_eit_destroy(&it);
    igraph_es_destroy(&es);
  }

  igraph_destroy(&g);

  return 0;
}
