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
#include <stdio.h>
#include <unistd.h>		/* unlink */

void custom_warning_handler (const char *reason, const char *file,
			     int line, int igraph_errno) {
  printf("Warning: %s\n", reason);
}

int main(int argc, char **argv) {
  igraph_t g;
  igraph_error_handler_t* oldhandler;
  igraph_warning_handler_t* oldwarnhandler;
  int result;
  FILE *ifile, *ofile;

  /* GraphML */
  ifile=fopen("test.gxl", "r");
  if (ifile==0) {
    return 10;
  }
  
  oldhandler=igraph_set_error_handler(igraph_error_handler_ignore);
  oldwarnhandler=igraph_set_warning_handler(custom_warning_handler);
  if ((result=igraph_read_graph_graphml(&g, ifile, 0))) {
    // maybe it is simply disabled at compile-time
    if (result == IGRAPH_UNIMPLEMENTED) return 77;
    return 1;
  }
  igraph_set_error_handler(oldhandler);

  fclose(ifile);

  /* Write it back */
  ofile=fopen("test2.gxl", "w");
  /* If we can't create the test file, just skip the test */
  if (ofile) {
    if ((result=igraph_write_graph_graphml(&g, ofile, /*prefixattr=*/ 1))) {
      return 1;
    }
    fclose(ofile);
    unlink("test2.gxl");
  }
  
  printf("The directed graph:\n");
  printf("Vertices: %li\n", (long int) igraph_vcount(&g));
  printf("Edges: %li\n", (long int) igraph_ecount(&g));
  printf("Directed: %i\n", (int) igraph_is_directed(&g));
  igraph_write_graph_edgelist(&g, stdout);
  igraph_destroy(&g);
 
  /* The same with undirected graph */
  ifile=fopen("test.gxl", "r");
  if ((result=igraph_read_graph_graphml(&g, ifile, 0))) {
    return 1;
  }
  fclose(ifile);
  printf("The undirected graph:\n");
  printf("Vertices: %li\n", (long int) igraph_vcount(&g));
  printf("Edges: %li\n", (long int) igraph_ecount(&g));
  printf("Directed: %i\n", (int) igraph_is_directed(&g));
  igraph_write_graph_edgelist(&g, stdout);
  igraph_destroy(&g);
  igraph_set_warning_handler(oldwarnhandler);

  /* There were sometimes problems with this file */
  /* Only if called from R though, and only on random occasions, once in every 
     ten reads. Do testing here doesn't make much sense, but if we have the file 
     then let's do it anyway. */
  ifile=fopen("graphml-hsa05010.xml", "r");  
  igraph_read_graph_graphml(&g, ifile, 0);
  fclose(ifile);
  igraph_destroy(&g);
  
  return 0;
}
