/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2021  The igraph development team

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

#if defined(VECTOR_LIST)
  /* It was indicated that every item in a list is a vector of the base type
   * so let's define ITEM_TYPE appropriately */
  #define ITEM_TYPE BASE_VECTOR
#elif defined(MATRIX_LIST)
  /* It was indicated that every item in a list is a matrix of the base type
   * so let's define ITEM_TYPE appropriately */
  #define ITEM_TYPE BASE_MATRIX
#else
  #define ITEM_TYPE BASE
#endif

/**
 * Vector list, dealing with lists of typed vectors efficiently.
 * \ingroup types
 */

typedef struct {
    ITEM_TYPE* stor_begin;
    ITEM_TYPE* stor_end;
    ITEM_TYPE* end;
#ifdef EXTRA_TYPE_FIELDS
    EXTRA_TYPE_FIELDS
#endif
} TYPE;

/*--------------------*/
/* Allocation         */
/*--------------------*/

IGRAPH_EXPORT igraph_error_t FUNCTION(init)(TYPE* v, igraph_integer_t size);
IGRAPH_EXPORT void FUNCTION(destroy)(TYPE* v);

/*--------------------*/
/* Accessing elements */
/*--------------------*/

IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE ITEM_TYPE* FUNCTION(get_ptr)(const TYPE* v, igraph_integer_t pos);
IGRAPH_EXPORT void FUNCTION(set)(TYPE* v, igraph_integer_t pos, ITEM_TYPE* e);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE ITEM_TYPE* FUNCTION(tail_ptr)(const TYPE *v);

/*-----------------*/
/* List properties */
/*-----------------*/

IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_integer_t FUNCTION(capacity)(const TYPE* v);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(empty)(const TYPE* v);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_integer_t FUNCTION(size)(const TYPE* v);

/*------------------------*/
/* Resizing operations    */
/*------------------------*/

IGRAPH_EXPORT void FUNCTION(clear)(TYPE* v);
IGRAPH_EXPORT igraph_error_t FUNCTION(reserve)(TYPE* v, igraph_integer_t capacity);
IGRAPH_EXPORT igraph_error_t FUNCTION(resize)(TYPE* v, igraph_integer_t new_size);

/*------------------------*/
/* Adding/removing items  */
/*------------------------*/

IGRAPH_EXPORT void FUNCTION(discard)(TYPE* v, igraph_integer_t index);
IGRAPH_EXPORT void FUNCTION(discard_back)(TYPE* v);
IGRAPH_EXPORT void FUNCTION(discard_fast)(TYPE* v, igraph_integer_t index);
IGRAPH_EXPORT igraph_error_t FUNCTION(insert)(TYPE* v, igraph_integer_t pos, ITEM_TYPE* e);
IGRAPH_EXPORT igraph_error_t FUNCTION(insert_copy)(TYPE* v, igraph_integer_t pos, const ITEM_TYPE* e);
IGRAPH_EXPORT igraph_error_t FUNCTION(insert_new)(TYPE* v, igraph_integer_t pos, ITEM_TYPE** result);
IGRAPH_EXPORT igraph_error_t FUNCTION(push_back)(TYPE* v, ITEM_TYPE* e);
IGRAPH_EXPORT igraph_error_t FUNCTION(push_back_copy)(TYPE* v, const ITEM_TYPE* e);
IGRAPH_EXPORT igraph_error_t FUNCTION(push_back_new)(TYPE* v, ITEM_TYPE** result);
IGRAPH_EXPORT ITEM_TYPE FUNCTION(pop_back)(TYPE* v);
IGRAPH_EXPORT igraph_error_t FUNCTION(remove)(TYPE* v, igraph_integer_t index, ITEM_TYPE* e);
IGRAPH_EXPORT igraph_error_t FUNCTION(remove_fast)(TYPE* v, igraph_integer_t index, ITEM_TYPE* e);
IGRAPH_EXPORT void FUNCTION(replace)(TYPE* v, igraph_integer_t pos, ITEM_TYPE* e);
IGRAPH_EXPORT void FUNCTION(remove_consecutive_duplicates)(TYPE *v, igraph_bool_t (*eq)(const ITEM_TYPE*, const ITEM_TYPE*));

/*------------------*/
/* Exchanging items */
/*------------------*/

IGRAPH_EXPORT igraph_error_t FUNCTION(permute)(TYPE *v, const igraph_vector_int_t *index);
IGRAPH_EXPORT igraph_error_t FUNCTION(reverse)(TYPE *v);
IGRAPH_EXPORT igraph_error_t FUNCTION(swap)(TYPE *v1, TYPE *v2);
IGRAPH_EXPORT igraph_error_t FUNCTION(swap_elements)(TYPE* v, igraph_integer_t i, igraph_integer_t j);

/*-----------*/
/* Sorting   */
/*-----------*/

IGRAPH_EXPORT void FUNCTION(sort)(
        TYPE *v, int (*cmp)(const ITEM_TYPE*, const ITEM_TYPE*));
IGRAPH_EXPORT igraph_error_t FUNCTION(sort_ind)(
        TYPE *v, igraph_vector_int_t *ind,
        int (*cmp)(const ITEM_TYPE*, const ITEM_TYPE*)
);

#undef ITEM_TYPE
