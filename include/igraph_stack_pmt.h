/*
   igraph library.
   Copyright (C) 2007-2025  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>

/**
 * Stack data type.
 * \ingroup internal
 */

typedef struct TYPE(igraph_stack) {
    BASE* stor_begin;
    BASE* stor_end;
    BASE* end;
} TYPE(igraph_stack);

IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_stack, init)(TYPE(igraph_stack)* s, igraph_int_t capacity);
IGRAPH_EXPORT void FUNCTION(igraph_stack, destroy)(TYPE(igraph_stack)* s);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_stack, reserve)(TYPE(igraph_stack)* s, igraph_int_t capacity);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_stack, empty)(TYPE(igraph_stack)* s);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_int_t FUNCTION(igraph_stack, size)(const TYPE(igraph_stack)* s);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_int_t FUNCTION(igraph_stack, capacity)(const TYPE(igraph_stack)* s);
IGRAPH_EXPORT void FUNCTION(igraph_stack, clear)(TYPE(igraph_stack)* s);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_stack, push)(TYPE(igraph_stack)* s, BASE elem);
IGRAPH_EXPORT BASE FUNCTION(igraph_stack, pop)(TYPE(igraph_stack)* s);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE BASE FUNCTION(igraph_stack, top)(const TYPE(igraph_stack)* s);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_stack, print)(const TYPE(igraph_stack)* s);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_stack, fprint)(const TYPE(igraph_stack)* s, FILE *file);
