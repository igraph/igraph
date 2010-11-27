/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2010  Gabor Csardi <csardi.gabor@gmail.com>
   Rue de l'Industrie 5, Lausanne 1005, Switzerland
   
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
  FILE *ifile;
  
  ifile=fopen("pajek5.net", "r");
  if (!ifile) { return 1; }
  igraph_read_graph_pajek(&g, ifile);
  fclose(ifile);
  if (igraph_vcount(&g) != 10 || igraph_ecount(&g) != 9 ||
      igraph_is_directed(&g)) { return 2; }
  igraph_destroy(&g);

  ifile=fopen("pajek6.net", "r");
  if (!ifile) { return 3; }
  igraph_read_graph_pajek(&g, ifile);
  fclose(ifile);
  if (igraph_vcount(&g) != 10 || igraph_ecount(&g) != 9 ||
      !igraph_is_directed(&g)) { return 3; }
  igraph_destroy(&g);
  
  return 0;
}
