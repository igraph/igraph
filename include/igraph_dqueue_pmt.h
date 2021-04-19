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

IGRAPH_EXPORT int FUNCTION(igraph_dqueue, init)    (TYPE(igraph_dqueue)* q, long int size);
IGRAPH_EXPORT void FUNCTION(igraph_dqueue, destroy) (TYPE(igraph_dqueue)* q);
IGRAPH_EXPORT igraph_bool_t FUNCTION(igraph_dqueue, empty)   (const TYPE(igraph_dqueue)* q);
IGRAPH_EXPORT void FUNCTION(igraph_dqueue, clear)   (TYPE(igraph_dqueue)* q);
IGRAPH_EXPORT igraph_bool_t FUNCTION(igraph_dqueue, full)    (TYPE(igraph_dqueue)* q);
IGRAPH_EXPORT long int FUNCTION(igraph_dqueue, size)    (const TYPE(igraph_dqueue)* q);
IGRAPH_EXPORT BASE FUNCTION(igraph_dqueue, pop)     (TYPE(igraph_dqueue)* q);
IGRAPH_EXPORT BASE FUNCTION(igraph_dqueue, pop_back)(TYPE(igraph_dqueue)* q);
IGRAPH_EXPORT BASE FUNCTION(igraph_dqueue, head)    (const TYPE(igraph_dqueue)* q);
IGRAPH_EXPORT BASE FUNCTION(igraph_dqueue, back)    (const TYPE(igraph_dqueue)* q);
IGRAPH_EXPORT int FUNCTION(igraph_dqueue, push)    (TYPE(igraph_dqueue)* q, BASE elem);
IGRAPH_EXPORT int FUNCTION(igraph_dqueue, print)(const TYPE(igraph_dqueue)* q);
IGRAPH_EXPORT int FUNCTION(igraph_dqueue, fprint)(const TYPE(igraph_dqueue)* q, FILE *file);
IGRAPH_EXPORT BASE FUNCTION(igraph_dqueue, e)(const TYPE(igraph_dqueue) *q, long int idx);
