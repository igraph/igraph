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

#include "igraph_types.h"
#include "igraph_matrix.h"

#define BASE_IGRAPH_REAL
#include "igraph_pmt.h"
#include "matrix.pmt"
#include "igraph_pmt_off.h"
#undef BASE_IGRAPH_REAL

#define BASE_LONG
#include "igraph_pmt.h"
#include "matrix.pmt"
#include "igraph_pmt_off.h"
#undef BASE_LONG

#define BASE_CHAR
#include "igraph_pmt.h"
#include "matrix.pmt"
#include "igraph_pmt_off.h"
#undef BASE_CHAR

#define BASE_BOOL
#include "igraph_pmt.h"
#include "matrix.pmt"
#include "igraph_pmt_off.h"
#undef BASE_BOOL

#define BASE_COMPLEX
#include "igraph_pmt.h"
#include "matrix.pmt"
#include "igraph_pmt_off.h"
#undef BASE_COMPLEX

int igraph_matrix_complex_print(const igraph_matrix_complex_t *m) {

  long int nr=igraph_matrix_complex_nrow(m);
  long int nc=igraph_matrix_complex_ncol(m);
  long int i, j;
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) {
      igraph_complex_t z=MATRIX(*m, i, j);
      if (j!=0) { putchar(' '); }
      if (IGRAPH_IMAG(z) < 0) {
	printf("%g-%gi", IGRAPH_REAL(z), IGRAPH_IMAG(z));
      } else {
	printf("%g+%gi", IGRAPH_REAL(z), IGRAPH_IMAG(z));
      }
    }
    printf("\n");
  }
  
  return 0;
}

int igraph_matrix_complex_fprint(const igraph_matrix_complex_t *m, 
				 FILE *file) {

  long int nr=igraph_matrix_complex_nrow(m);
  long int nc=igraph_matrix_complex_ncol(m);
  long int i, j;
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) {
      igraph_complex_t z=MATRIX(*m, i, j);
      if (j!=0) { putchar(' '); }
      if (IGRAPH_IMAG(z) < 0) {
	fprintf(file, "%g-%gi", IGRAPH_REAL(z), IGRAPH_IMAG(z));
      } else {
	fprintf(file, "%g+%gi", IGRAPH_REAL(z), IGRAPH_IMAG(z));
      }
    }
    fprintf(file, "\n");
  }
  
  return 0;
}
