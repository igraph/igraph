/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2014  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge, MA 02139, USA

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

  igraph_vector_t longvec;
  igraph_vector_t idx, idx2;

#define N 10

  igraph_vector_init(&longvec, N);
  igraph_vector_init_seq(&idx, 0, N-1);
  igraph_vector_init(&idx2, 0);

  igraph_vector_qsort_ind_stable(&longvec, &idx2, /* decreasing= */ 0);
  if (!igraph_vector_all_e(&idx, &idx2)) { return 1; }

  igraph_vector_null(&idx2);
  igraph_vector_clear(&idx2);
  igraph_vector_qsort_ind_stable(&longvec, &idx2, /* decreasing= */ 1);
  if (!igraph_vector_all_e(&idx, &idx2)) { return 2; }

  igraph_vector_destroy(&idx);
  igraph_vector_destroy(&idx2);
  igraph_vector_destroy(&longvec);

  return 0;
}
