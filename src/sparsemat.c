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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include "sparsemat.h"
#include "error.h"

int igraph_sparsemat_init(igraph_sparsemat_t *A, int rows, int cols, int nzmax) {

  if (rows < 0) { 
    IGRAPH_ERROR("Negative number of rows", IGRAPH_EINVAL);
  }
  if (cols < 0) {
    IGRAPH_ERROR("Negative number of columns", IGRAPH_EINVAL);
  }
  
  A->cs=cs_spalloc( rows, cols, nzmax, /*values=*/ 1, 
		  /*triplet=*/ 1);
  if (!A->cs) {
    IGRAPH_ERROR("Cannot allocate memory for sparse matrix", IGRAPH_ENOMEM);
  }

  return 0;
}

void igraph_sparsemat_destroy(igraph_sparsemat_t *A) {
  cs_spfree(A->cs);
}

int igraph_sparsemat_realloc(igraph_sparsemat_t *A, int nzmax) {
  return !cs_sprealloc(A->cs, nzmax);
}
