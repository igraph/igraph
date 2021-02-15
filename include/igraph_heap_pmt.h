/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

typedef struct TYPE(igraph_heap) {
    BASE* stor_begin;
    BASE* stor_end;
    BASE* end;
    int destroy;
} TYPE(igraph_heap);

IGRAPH_EXPORT int FUNCTION(igraph_heap, init)(TYPE(igraph_heap)* h, long int size);
IGRAPH_EXPORT int FUNCTION(igraph_heap, init_array)(TYPE(igraph_heap) *t, BASE* data, long int len);
IGRAPH_EXPORT void FUNCTION(igraph_heap, destroy)(TYPE(igraph_heap)* h);
IGRAPH_EXPORT igraph_bool_t FUNCTION(igraph_heap, empty)(TYPE(igraph_heap)* h);
IGRAPH_EXPORT int FUNCTION(igraph_heap, push)(TYPE(igraph_heap)* h, BASE elem);
IGRAPH_EXPORT BASE FUNCTION(igraph_heap, top)(TYPE(igraph_heap)* h);
IGRAPH_EXPORT BASE FUNCTION(igraph_heap, delete_top)(TYPE(igraph_heap)* h);
IGRAPH_EXPORT long int FUNCTION(igraph_heap, size)(TYPE(igraph_heap)* h);
IGRAPH_EXPORT int FUNCTION(igraph_heap, reserve)(TYPE(igraph_heap)* h, long int size);
