/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

int main() {

  igraph_t g1, g2, res;
  igraph_vector_t v;

  /* composition with the empty graph */
  igraph_empty(&g1, 5, IGRAPH_DIRECTED);
  igraph_full(&g2, 5, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
  igraph_compose(&res, &g1, &g2);
  if (igraph_ecount(&res) != 0) { 
    return 1;
  }
  igraph_destroy(&res);
  igraph_compose(&res, &g2, &g1);
  if (igraph_ecount(&res) != 0) { 
    return 2;
  }
  igraph_destroy(&res);
  igraph_destroy(&g1);
  igraph_destroy(&g2);

  /* same but undirected */
  igraph_empty(&g1, 5, IGRAPH_UNDIRECTED);
  igraph_full(&g2, 5, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
  igraph_compose(&res, &g1, &g2);
  if (igraph_ecount(&res) != 0) { 
    return 1;
  }
  igraph_destroy(&res);
  igraph_compose(&res, &g2, &g1);
  if (igraph_ecount(&res) != 0) { 
    return 2;
  }
  igraph_destroy(&res);
  igraph_destroy(&g1);
  igraph_destroy(&g2);

  /* proper directed graph */
  igraph_vector_init_int_end(&v, -1, 0,1, 1,2, 5,6, -1);
  igraph_create(&g1, &v, 0, IGRAPH_DIRECTED);
  igraph_vector_destroy(&v);

  igraph_vector_init_int_end(&v, -1, 0,1, 2,4, 5,6, -1);
  igraph_create(&g2, &v, 0, IGRAPH_DIRECTED);
  igraph_vector_destroy(&v);
  
  igraph_compose(&res, &g1, &g2);
  igraph_write_graph_edgelist(&res, stdout);
  igraph_destroy(&res);
  igraph_destroy(&g1);
  igraph_destroy(&g2);

  /* undirected graph */
  igraph_vector_init_int_end(&v, -1, 0,1, 1,2, 5,6, -1);
  igraph_create(&g1, &v, 0, IGRAPH_UNDIRECTED);
  igraph_vector_destroy(&v);

  igraph_vector_init_int_end(&v, -1, 0,1, 0,4, 5,6, -1);
  igraph_create(&g2, &v, 0, IGRAPH_UNDIRECTED);
  igraph_vector_destroy(&v);
  
  igraph_compose(&res, &g1, &g2);
  igraph_write_graph_edgelist(&res, stdout);
  igraph_destroy(&res);
  igraph_destroy(&g1);
  igraph_destroy(&g2);

  return 0;
}
