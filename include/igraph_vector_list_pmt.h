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

/**
 * Vector list, dealing with lists of typed vectors efficiently.
 * \ingroup types
 */

typedef struct {
    BASE_VECTOR* stor_begin;
    BASE_VECTOR* stor_end;
    BASE_VECTOR* end;
} TYPE(igraph_vector);

/*--------------------*/
/* Allocation         */
/*--------------------*/

IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, init)(
        TYPE(igraph_vector)* v, igraph_integer_t size);
IGRAPH_EXPORT void FUNCTION(igraph_vector, destroy)(TYPE(igraph_vector)* v);

IGRAPH_EXPORT igraph_integer_t FUNCTION(igraph_vector, capacity)(const TYPE(igraph_vector)*v);

/*--------------------*/
/* Accessing elements */
/*--------------------*/

IGRAPH_EXPORT BASE_VECTOR* FUNCTION(igraph_vector, get)(const TYPE(igraph_vector)* v, igraph_integer_t pos);
IGRAPH_EXPORT void FUNCTION(igraph_vector, replace)(TYPE(igraph_vector)* v, igraph_integer_t pos, BASE_VECTOR* value);
IGRAPH_EXPORT void FUNCTION(igraph_vector, set)(TYPE(igraph_vector)* v, igraph_integer_t pos, BASE_VECTOR* value);
IGRAPH_EXPORT BASE_VECTOR* FUNCTION(igraph_vector, tail)(const TYPE(igraph_vector) *v);

/*-----------------*/
/* List properties */
/*-----------------*/

IGRAPH_EXPORT igraph_bool_t FUNCTION(igraph_vector, empty)(const TYPE(igraph_vector)* v);
IGRAPH_EXPORT igraph_integer_t FUNCTION(igraph_vector, size)(const TYPE(igraph_vector)* v);

/*------------------------*/
/* Resizing operations    */
/*------------------------*/

IGRAPH_EXPORT void FUNCTION(igraph_vector, clear)(TYPE(igraph_vector)* v);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, resize)(
        TYPE(igraph_vector)* v, igraph_integer_t new_size);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, reserve)(
        TYPE(igraph_vector)* v, igraph_integer_t capacity);

/*------------------------*/
/* Adding/removing items  */
/*------------------------*/

IGRAPH_EXPORT void FUNCTION(igraph_vector, discard)(
        TYPE(igraph_vector)* v, igraph_integer_t index);
IGRAPH_EXPORT void FUNCTION(igraph_vector, discard_back)(TYPE(igraph_vector)* v);
IGRAPH_EXPORT void FUNCTION(igraph_vector, discard_fast)(
        TYPE(igraph_vector)* v, igraph_integer_t index);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, push_back)(
        TYPE(igraph_vector)* v, BASE_VECTOR* e);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, push_back_copy)(
        TYPE(igraph_vector)* v, BASE_VECTOR* e);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, push_back_new)(
        TYPE(igraph_vector)* v, BASE_VECTOR** result);
IGRAPH_EXPORT BASE_VECTOR* FUNCTION(igraph_vector, pop_back)(TYPE(igraph_vector)* v);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, remove)(
        TYPE(igraph_vector)* v, igraph_integer_t index, BASE_VECTOR* result);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, remove_fast)(
        TYPE(igraph_vector)* v, igraph_integer_t index, BASE_VECTOR* result);

/*-----------*/
/* Sorting   */
/*-----------*/

IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, permute)(
        TYPE(igraph_vector) *v, const igraph_vector_int_t *index);
IGRAPH_EXPORT void FUNCTION(igraph_vector, sort)(
        TYPE(igraph_vector) *v, int (*cmp)(const BASE_VECTOR*, const BASE_VECTOR*));
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_vector, sort_ind)(
        TYPE(igraph_vector) *v, igraph_vector_int_t *ind,
        int (*cmp)(const BASE_VECTOR*, const BASE_VECTOR*)
);
