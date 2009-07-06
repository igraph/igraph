/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2009  Gabor Csardi <Gabor.Csardi@unil.ch>
   UNIL DGM, Rue de Bugnon 27, CH-1005 Lausanne, Switzerland
   
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include <igraph.h>
#include <stdio.h>

int main() {

  igraph_t g;
  FILE *karate, *neural;
  igraph_real_t res;
  igraph_vector_long_t types;
  igraph_vector_t degree, outdegree, indegree;
  
  karate=fopen("karate.gml", "r");
  igraph_read_graph_gml(&g, karate);
  fclose(karate);
  
  igraph_vector_long_init(&types, 0);
  igraph_vector_init(&degree, 0);
  igraph_degree(&g, &degree, igraph_vss_all(), IGRAPH_ALL, /*loops=*/ 1);
  igraph_vector_floor(&degree, &types);

  igraph_assortativity_nominal(&g, &types, &res, /*directed=*/ 0);
  printf("%.5f\n", res);
  
  igraph_destroy(&g);

  /*---------------------*/

  neural=fopen("celegansneural.gml", "r");
  igraph_read_graph_gml(&g, neural);
  fclose(neural);
  
  igraph_degree(&g, &degree, igraph_vss_all(), IGRAPH_ALL, /*loops=*/ 1);
  igraph_vector_floor(&degree, &types);

  igraph_assortativity_nominal(&g, &types, &res, /*directed=*/ 1);
  printf("%.5f\n", res);  
  igraph_assortativity_nominal(&g, &types, &res, /*directed=*/ 0);
  printf("%.5f\n", res);  

  igraph_destroy(&g);
  igraph_vector_long_destroy(&types);
  igraph_vector_destroy(&degree);

  /*---------------------*/
  
  karate=fopen("karate.gml", "r");
  igraph_read_graph_gml(&g, karate);
  fclose(karate);
  
  igraph_vector_init(&degree, 0);
  igraph_degree(&g, &degree, igraph_vss_all(), IGRAPH_ALL, /*loops=*/ 1);
  igraph_vector_add_constant(&degree, -1);

  igraph_assortativity(&g, &degree, 0, &res, /*directed=*/ 0);
  printf("%.5f\n", res);
  
  igraph_destroy(&g);

  /*---------------------*/

  neural=fopen("celegansneural.gml", "r");
  igraph_read_graph_gml(&g, neural);
  fclose(neural);
  
  igraph_degree(&g, &degree, igraph_vss_all(), IGRAPH_ALL, /*loops=*/ 1);
  igraph_vector_add_constant(&degree, -1);

  igraph_assortativity(&g, &degree, 0, &res, /*directed=*/ 1);
  printf("%.5f\n", res);  
  igraph_assortativity(&g, &degree, 0, &res, /*directed=*/ 0);
  printf("%.5f\n", res);  

  igraph_vector_destroy(&degree);

  /*---------------------*/
  
  igraph_vector_init(&indegree, 0);
  igraph_vector_init(&outdegree, 0);
  igraph_degree(&g, &indegree, igraph_vss_all(), IGRAPH_IN, /*loops=*/ 1);
  igraph_degree(&g, &outdegree, igraph_vss_all(), IGRAPH_OUT, /*loops=*/ 1);
  igraph_vector_add_constant(&indegree, -1);
  igraph_vector_add_constant(&outdegree, -1);
  
  igraph_assortativity(&g, &outdegree, &indegree, &res, /*directed=*/ 1);
  printf("%.5f\n", res);

  igraph_vector_destroy(&indegree);
  igraph_vector_destroy(&outdegree);

  /*---------------------*/
  
  igraph_assortativity_degree(&g, &res, /*directed=*/ 1);
  printf("%.5f\n", res);

  igraph_destroy(&g);  

  /*---------------------*/

  karate=fopen("karate.gml", "r");
  igraph_read_graph_gml(&g, karate);
  fclose(karate);
  
  igraph_assortativity_degree(&g, &res, /*directed=*/ 1);
  printf("%.5f\n", res);
  
  igraph_destroy(&g);
  
  return 0;
}

