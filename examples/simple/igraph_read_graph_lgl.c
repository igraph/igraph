/* -*- mode: C -*-  */
/* 
   IGraph R package.
   Copyright (C) 2005 Gabor Csardi <csardi@rmki.kfki.hu>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include <igraph.h>

int main() {

  igraph_t g;
  FILE *input;

  /* Without names and weights */
  input=fopen("igraph_read_graph_lgl-1.lgl", "r");
  if (!input) { 
    return 1;
  }
  igraph_read_graph_lgl(&g, input, 0, 0);
  fclose(input);
  igraph_write_graph_edgelist(&g, stdout);
  igraph_destroy(&g);

  /* With names and weights */
  input=fopen("igraph_read_graph_lgl-2.lgl", "r");
  if (!input) {
    return 2;
  }
  igraph_read_graph_lgl(&g, input, 0, 0);
  fclose(input);
  igraph_write_graph_ncol(&g, stdout, 0, 0);
  igraph_destroy(&g);

  /* Erroneous LGL file (empty vertex name) */
  input=fopen("igraph_read_graph_lgl-3.lgl", "r");
  if (!input) {
    return 3;
  }
  igraph_set_error_handler(igraph_error_handler_ignore);
  if (igraph_read_graph_lgl(&g, input, 0, 0) != IGRAPH_PARSEERROR) {
    return 4;
  }
  fclose(input);

  return 0;
}
