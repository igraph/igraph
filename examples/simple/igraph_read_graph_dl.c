/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2009  Gabor Csardi <csardi@rmki.kfki.hu>
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
#include <stdio.h>
#include <stdlib.h>

int main() {
  
  const char *files[] = { "fullmatrix1.dl", "fullmatrix2.dl", 
			  "fullmatrix3.dl", "fullmatrix4.dl",
			  "edgelist1.dl", "edgelist2.dl", "edgelist3.dl",
			  "nodelist1.dl", "nodelist2.dl" };
  int no_files=sizeof(files)/sizeof(const char*);
  int i, ret;
  igraph_t g;
  FILE *infile;

  for (i=0; i<no_files; i++) {
    infile=fopen(files[i], "r");
    if (!infile) {
      printf("Cannot open file: %s\n", files[i]);
      exit(1+i);
    }
    igraph_read_graph_dl(&g, infile);
    ret=fclose(infile);
    if (ret) {
      printf("Cannot close file: %s\n", files[i]);
      exit(11+i);
    }
/*     igraph_destroy(&g); */
  }

  return 0;
}
