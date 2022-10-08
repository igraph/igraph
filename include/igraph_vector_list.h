/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2022  The igraph development team

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

#ifndef IGRAPH_VECTOR_LIST_H
#define IGRAPH_VECTOR_LIST_H

#include "igraph_decls.h"
#include "igraph_types.h"
#include "igraph_vector.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Flexible list of vectors                           */
/* -------------------------------------------------- */

/* Indicate to igraph_typed_list_pmt.h that we are going to work with _vectors_
 * of the base type, not the base type directly */
#define VECTOR_LIST

#define BASE_IGRAPH_REAL
#include "igraph_pmt.h"
#include "igraph_typed_list_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_IGRAPH_REAL

#define BASE_INT
#include "igraph_pmt.h"
#include "igraph_typed_list_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_INT

#undef VECTOR_LIST

/* -------------------------------------------------- */
/* Helper macros                                      */
/* -------------------------------------------------- */

#define IGRAPH_VECTOR_LIST_INIT_FINALLY(v, size) \
    do { IGRAPH_CHECK(igraph_vector_list_init(v, size)); \
        IGRAPH_FINALLY(igraph_vector_list_destroy, v); } while (0)
#define IGRAPH_VECTOR_BOOL_LIST_INIT_FINALLY(v, size) \
    do { IGRAPH_CHECK(igraph_vector_bool_list_init(v, size)); \
        IGRAPH_FINALLY(igraph_vector_bool_list_destroy, v); } while (0)
#define IGRAPH_VECTOR_CHAR_LIST_INIT_FINALLY(v, size) \
  do { IGRAPH_CHECK(igraph_vector_char_list_init(v, size)); \
  IGRAPH_FINALLY(igraph_vector_char_list_destroy, v); } while (0)
#define IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(v, size) \
    do { IGRAPH_CHECK(igraph_vector_int_list_init(v, size)); \
        IGRAPH_FINALLY(igraph_vector_int_list_destroy, v); } while (0)

__END_DECLS

#endif
