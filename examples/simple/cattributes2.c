/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2007  Gabor Csardi <csardi@rmki.kfki.hu>
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
  igraph_vector_t y;  

  /* turn on attribute handling */
  igraph_i_set_attribute_table(&igraph_cattribute_table);

  /* Create a graph, add some attributes and save it as a GraphML file */
  igraph_famous(&g, "Petersen");
  SETGAS(&g, "name", "Petersen's graph");
  SETGAN(&g, "vertices", igraph_vcount(&g));
  SETGAN(&g, "edges", igraph_ecount(&g));

  igraph_vector_init_seq(&y, 1, igraph_vcount(&g));
  SETVANV(&g, "id", &y);
  igraph_vector_destroy(&y);

  SETVAS(&g, "name", 0, "foo");
  SETVAS(&g, "name", 1, "foobar");
  
  igraph_vector_init_seq(&y, 1, igraph_ecount(&g));
  SETEANV(&g, "id", &y);
  igraph_vector_destroy(&y);

  SETEAS(&g, "name", 0, "FOO");
  SETEAS(&g, "name", 1, "FOOBAR");

  igraph_write_graph_gml(&g, stdout, 0, "");
  igraph_write_graph_graphml(&g, stdout);
   
  igraph_destroy(&g);
  
  return 0;
}
