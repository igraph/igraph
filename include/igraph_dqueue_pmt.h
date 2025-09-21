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

/**
 * Double ended queue data type.
 * \ingroup internal
 */

typedef struct TYPE(igraph_dqueue) {
    BASE *begin;
    BASE *end;
    BASE *stor_begin;
    BASE *stor_end;
} TYPE(igraph_dqueue);

IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_dqueue, init)(TYPE(igraph_dqueue)* q, igraph_int_t capacity);
IGRAPH_EXPORT void FUNCTION(igraph_dqueue, destroy)(TYPE(igraph_dqueue)* q);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_dqueue, empty)(const TYPE(igraph_dqueue)* q);
IGRAPH_EXPORT void FUNCTION(igraph_dqueue, clear)(TYPE(igraph_dqueue)* q);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_dqueue, full)(TYPE(igraph_dqueue)* q);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_int_t FUNCTION(igraph_dqueue, size)(const TYPE(igraph_dqueue)* q);
IGRAPH_EXPORT BASE FUNCTION(igraph_dqueue, pop)(TYPE(igraph_dqueue)* q);
IGRAPH_EXPORT BASE FUNCTION(igraph_dqueue, pop_back)(TYPE(igraph_dqueue)* q);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE BASE FUNCTION(igraph_dqueue, head)(const TYPE(igraph_dqueue)* q);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE BASE FUNCTION(igraph_dqueue, back)(const TYPE(igraph_dqueue)* q);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_dqueue, push)(TYPE(igraph_dqueue)* q, BASE elem);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_dqueue, print)(const TYPE(igraph_dqueue)* q);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_dqueue, fprint)(const TYPE(igraph_dqueue)* q, FILE *file);
IGRAPH_EXPORT BASE FUNCTION(igraph_dqueue, get)(const TYPE(igraph_dqueue) *q, igraph_int_t idx);
