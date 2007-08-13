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

#include "types.h"

#define BASE_IGRAPH_REAL
#include "igraph_pmt.h"
#include "vector.pmt"
#include "igraph_pmt_off.h"
#undef BASE_IGRAPH_REAL

#define BASE_LONG
#include "igraph_pmt.h"
#include "vector.pmt"
#include "igraph_pmt_off.h"
#undef BASE_LONG

#define BASE_CHAR
#include "igraph_pmt.h"
#include "vector.pmt"
#include "igraph_pmt_off.h"
#undef BASE_CHAR

#define BASE_BOOL
#include "igraph_pmt.h"
#include "vector.pmt"
#include "igraph_pmt_off.h"
#undef BASE_BOOL

int igraph_vector_order2(igraph_vector_t *v) {

  igraph_indheap_t heap;
  
  igraph_indheap_init_array(&heap, VECTOR(*v), igraph_vector_size(v));
  IGRAPH_FINALLY(igraph_indheap_destroy, &heap);

  igraph_vector_clear(v);
  while (!igraph_indheap_empty(&heap)) {
    IGRAPH_CHECK(igraph_vector_push_back(v, igraph_indheap_max_index(&heap)-1));
    igraph_indheap_delete_max(&heap);
  }
  
  igraph_indheap_destroy(&heap);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

