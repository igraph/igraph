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

int main() {

  IGRAPH_VECTOR_CONSTANT(v1, 1, 2, 3, 4, 5);
  IGRAPH_VECTOR_CONSTANT(v2, 1, 2, 3, 4, 5);
  IGRAPH_VECTOR_CONSTANT(v3);
  IGRAPH_VECTOR_CONSTANT(v4, -1.5);

  igraph_vector_print(&v1);
  igraph_vector_print(&v2);
  igraph_vector_print(&v3);
  igraph_vector_print(&v4);

  /* ----------------------- */

  IGRAPH_VECTOR_FLOAT_CONSTANT(f1, 1, 2, 3, 4, 5);
  IGRAPH_VECTOR_FLOAT_CONSTANT(f2, 1, 2, 3, 4, 5);
  IGRAPH_VECTOR_FLOAT_CONSTANT(f3);
  IGRAPH_VECTOR_FLOAT_CONSTANT(f4, -1.5);

  igraph_vector_float_print(&f1);
  igraph_vector_float_print(&f2);
  igraph_vector_float_print(&f3);
  igraph_vector_float_print(&f4);

  /* ----------------------- */

  IGRAPH_VECTOR_LONG_CONSTANT(l1, 1, 2, 3, 4, 5);
  IGRAPH_VECTOR_LONG_CONSTANT(l2, 1, 2, 3, 4, 5);
  IGRAPH_VECTOR_LONG_CONSTANT(l3);
  IGRAPH_VECTOR_LONG_CONSTANT(l4, -1);

  igraph_vector_long_print(&l1);
  igraph_vector_long_print(&l2);
  igraph_vector_long_print(&l3);
  igraph_vector_long_print(&l4);

  /* ----------------------- */

  IGRAPH_VECTOR_CHAR_CONSTANT(c1, 1, 2, 3, 4, 5);
  IGRAPH_VECTOR_CHAR_CONSTANT(c2, 1, 2, 3, 4, 5);
  IGRAPH_VECTOR_CHAR_CONSTANT(c3);
  IGRAPH_VECTOR_CHAR_CONSTANT(c4, -1);

  igraph_vector_char_print(&c1);
  igraph_vector_char_print(&c2);
  igraph_vector_char_print(&c3);
  igraph_vector_char_print(&c4);
  
  /* ----------------------- */

  IGRAPH_VECTOR_BOOL_CONSTANT(b1, 1, 0, 1, 0, 1);
  IGRAPH_VECTOR_BOOL_CONSTANT(b2, 0, 1, 0, 1, 0);
  IGRAPH_VECTOR_BOOL_CONSTANT(b3);
  IGRAPH_VECTOR_BOOL_CONSTANT(b4, 1);

  igraph_vector_bool_print(&b1);
  igraph_vector_bool_print(&b2);
  igraph_vector_bool_print(&b3);
  igraph_vector_bool_print(&b4);

  /* ----------------------- */

  IGRAPH_VECTOR_INT_CONSTANT(i1, 1, 2, 3, 4, 5);
  IGRAPH_VECTOR_INT_CONSTANT(i2, 1, 2, 3, 4, 5);
  IGRAPH_VECTOR_INT_CONSTANT(i3);
  IGRAPH_VECTOR_INT_CONSTANT(i4, -1);

  igraph_vector_int_print(&i1);
  igraph_vector_int_print(&i2);
  igraph_vector_int_print(&i3);
  igraph_vector_int_print(&i4);

  /* ----------------------- */

  IGRAPH_VECTOR_TIME_CONSTANT(t1, 1, 2, 3, 4, 5);
  IGRAPH_VECTOR_TIME_CONSTANT(t2, 1, 2, 3, 4, 5);
  IGRAPH_VECTOR_TIME_CONSTANT(t3);
  IGRAPH_VECTOR_TIME_CONSTANT(t4, 1);

  igraph_vector_time_print(&t1);
  igraph_vector_time_print(&t2);
  igraph_vector_time_print(&t3);
  igraph_vector_time_print(&t4);

  return 0;
}
