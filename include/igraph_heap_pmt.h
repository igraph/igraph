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

typedef struct TYPE(igraph_heap) {
    BASE* stor_begin;
    BASE* stor_end;
    BASE* end;
} TYPE(igraph_heap);

IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_heap, init)(TYPE(igraph_heap)* h, igraph_int_t capacity);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_heap, init_array)(TYPE(igraph_heap) *t, const BASE *data, igraph_int_t len);
IGRAPH_EXPORT void FUNCTION(igraph_heap, destroy)(TYPE(igraph_heap)* h);
IGRAPH_EXPORT void FUNCTION(igraph_heap, clear)(TYPE(igraph_heap)* h);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_heap, empty)(const TYPE(igraph_heap)* h);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_heap, push)(TYPE(igraph_heap)* h, BASE elem);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE BASE FUNCTION(igraph_heap, top)(const TYPE(igraph_heap)* h);
IGRAPH_EXPORT BASE FUNCTION(igraph_heap, delete_top)(TYPE(igraph_heap)* h);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_int_t FUNCTION(igraph_heap, size)(const TYPE(igraph_heap)* h);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_heap, reserve)(TYPE(igraph_heap)* h, igraph_int_t capacity);
