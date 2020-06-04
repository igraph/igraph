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

DECLDIR int FUNCTION(igraph_heap, init)(TYPE(igraph_heap)* h, long int size);
DECLDIR int FUNCTION(igraph_heap, init_array)(TYPE(igraph_heap) *t, BASE* data, long int len);
DECLDIR void FUNCTION(igraph_heap, destroy)(TYPE(igraph_heap)* h);
DECLDIR igraph_bool_t FUNCTION(igraph_heap, empty)(TYPE(igraph_heap)* h);
DECLDIR int FUNCTION(igraph_heap, push)(TYPE(igraph_heap)* h, BASE elem);
DECLDIR BASE FUNCTION(igraph_heap, top)(TYPE(igraph_heap)* h);
DECLDIR BASE FUNCTION(igraph_heap, delete_top)(TYPE(igraph_heap)* h);
DECLDIR long int FUNCTION(igraph_heap, size)(TYPE(igraph_heap)* h);
DECLDIR int FUNCTION(igraph_heap, reserve)(TYPE(igraph_heap)* h, long int size);

void FUNCTION(igraph_heap, i_build)(BASE* arr, long int size, long int head);
void FUNCTION(igraph_heap, i_shift_up)(BASE* arr, long int size, long int elem);
void FUNCTION(igraph_heap, i_sink)(BASE* arr, long int size, long int head);
void FUNCTION(igraph_heap, i_switch)(BASE* arr, long int e1, long int e2);

