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

DECLDIR int FUNCTION(igraph_stack, init)(TYPE(igraph_stack)* s, long int size);
DECLDIR void FUNCTION(igraph_stack, destroy)(TYPE(igraph_stack)* s);
DECLDIR int FUNCTION(igraph_stack, reserve)(TYPE(igraph_stack)* s, long int size);
DECLDIR igraph_bool_t FUNCTION(igraph_stack, empty)(TYPE(igraph_stack)* s);
DECLDIR long int FUNCTION(igraph_stack, size)(const TYPE(igraph_stack)* s);
DECLDIR void FUNCTION(igraph_stack, clear)(TYPE(igraph_stack)* s);
DECLDIR int FUNCTION(igraph_stack, push)(TYPE(igraph_stack)* s, BASE elem);
DECLDIR BASE FUNCTION(igraph_stack, pop)(TYPE(igraph_stack)* s);
DECLDIR BASE FUNCTION(igraph_stack, top)(const TYPE(igraph_stack)* s);
DECLDIR int FUNCTION(igraph_stack, print)(const TYPE(igraph_stack)* s);
DECLDIR int FUNCTION(igraph_stack, fprint)(const TYPE(igraph_stack)* s, FILE *file);
