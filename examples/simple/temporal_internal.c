/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2014  Gabor Csardi <csardi.gabor@gmail.com>
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
#include <stdlib.h>

int test_igraph_i_index_vector() {

  igraph_vector_time_t v;
  igraph_vector_int_t idx, idx_correct;
  const int max_elem = 5;

  igraph_vector_time_init_int_end(&v, /* endmark= */ -1,
				  0,0, 1, 3,3,3, 4,4,4,
				  -1);
  igraph_vector_int_init_int_end(&idx_correct, /* endmark= */ -1,
				 0, 2, 3, 3, 6, 9, 9,
				 -1);
  igraph_vector_int_init(&idx, 0);
  igraph_vector_time_i_index(&v, &idx, max_elem);
  if (!igraph_vector_int_all_e(&idx, &idx_correct)) {
    igraph_vector_int_print(&idx);
    exit(1);
  }

  igraph_vector_time_destroy(&v);
  igraph_vector_int_destroy(&idx_correct);
  igraph_vector_int_destroy(&idx);

  igraph_vector_time_init_int_end(&v, /*endmark= */ -1,
				  2,2,2, 4, 5,5,5,5,5,
				  -1);
  igraph_vector_int_init_int_end(&idx_correct, /* endmark= */ -1,
				 0, 0, 0, 3, 3, 4, 9,
				 -1);
  igraph_vector_int_init(&idx, 0);
  igraph_vector_time_i_index(&v, &idx, max_elem);
  if (!igraph_vector_int_all_e(&idx, &idx_correct)) {
    igraph_vector_int_print(&idx);
    exit(2);
  }

  igraph_vector_time_destroy(&v);
  igraph_vector_int_destroy(&idx_correct);
  igraph_vector_int_destroy(&idx);

  return 0;
}

int test_igraph_i_index_vector_2() {

  igraph_vector_int_t idx;
  igraph_vector_int_init(&idx, 0);
  IGRAPH_VECTOR_TIME_CONSTANT(birth, 0, 0, 0, 0, 0);
  IGRAPH_VECTOR_INT_CONSTANT(correct, 0, 5);

  igraph_vector_time_i_index(&birth, &idx, 0);
  if (!igraph_vector_int_all_e(&idx, &correct)) { exit(4); }
  igraph_vector_int_destroy(&idx);

  return 0;
}

int test_igraph_i_index_through() {

  igraph_vector_t v, th, idx, idx_correct;
  const double max_elem = 5.0;

  igraph_vector_init_int_end(&v, /* endmark= */ -1,
			     0, 3, 4, 3, 3, 1, 4, 4, 0,
			     -1);
  igraph_vector_init_int_end(&th, /*endmark= */ -1,
			     8, 0, 5, 4, 3, 1, 2, 6, 7,
			     -1);
  igraph_vector_init_int_end(&idx_correct, /* endmark=*/ -1,
			     0, 2, 3, 3, 6, 9, 9,
			     -1);
  igraph_vector_init(&idx, 0);
  igraph_vector_i_index_through(&v, &th, &idx, max_elem);
  if (!igraph_vector_all_e(&idx, &idx_correct)) {
    igraph_vector_print(&idx);
    exit(3);
  }

  igraph_vector_destroy(&v);
  igraph_vector_destroy(&idx_correct);
  igraph_vector_destroy(&idx);
  igraph_vector_destroy(&th);

  igraph_vector_init_int_end(&v, /* endmark= */ -1,
			     5, 5, 2, 2, 5, 5, 4, 5, 2,
			     -1);
  igraph_vector_init_int_end(&th, /*endmark= */ -1,
			     2, 3, 8, 6, 7, 4, 5, 0, 1,
			     -1);
  igraph_vector_init_int_end(&idx_correct, /* endmark=*/ -1,
			     0, 0, 0, 3, 3, 4, 9,
			     -1);
  igraph_vector_init(&idx, 0);
  igraph_vector_i_index_through(&v, &th, &idx, max_elem);
  if (!igraph_vector_all_e(&idx, &idx_correct)) {
    igraph_vector_print(&idx);
    exit(3);
  }

  igraph_vector_destroy(&v);
  igraph_vector_destroy(&idx_correct);
  igraph_vector_destroy(&idx);
  igraph_vector_destroy(&th);

  return 0;
}

int test_igraph_add_vertices_at() {

  igraph_t graph;

  /* TODO */

  return 0;
}

int main() {
  test_igraph_i_index_vector();
  test_igraph_i_index_vector_2();
  test_igraph_i_index_through();
  test_igraph_add_vertices_at();
  return 0;
}
