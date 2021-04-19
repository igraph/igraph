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

typedef struct TYPE(igraph_array3) {
    TYPE(igraph_vector) data;
    long int n1, n2, n3, n1n2;
} TYPE(igraph_array3);

#ifndef IGRAPH_ARRAY3_INIT_FINALLY
#define IGRAPH_ARRAY3_INIT_FINALLY(a, n1, n2, n3) \
    do { IGRAPH_CHECK(igraph_array3_init(a, n1, n2, n3)); \
        IGRAPH_FINALLY(igraph_array3_destroy, a); } while (0)
#endif

#ifndef ARRAY3
    #define ARRAY3(m,i,j,k) ((m).data.stor_begin[(m).n1n2*(k)+(m).n1*(j)+(i)])
#endif

IGRAPH_EXPORT int FUNCTION(igraph_array3, init)(TYPE(igraph_array3) *a, long int n1, long int n2,
                                                long int n3);
IGRAPH_EXPORT void FUNCTION(igraph_array3, destroy)(TYPE(igraph_array3) *a);
IGRAPH_EXPORT long int FUNCTION(igraph_array3, size)(const TYPE(igraph_array3) *a);
IGRAPH_EXPORT long int FUNCTION(igraph_array3, n)(const TYPE(igraph_array3) *a, long int idx);
IGRAPH_EXPORT int FUNCTION(igraph_array3, resize)(TYPE(igraph_array3) *a, long int n1, long int n2,
                                                  long int n3);
IGRAPH_EXPORT void FUNCTION(igraph_array3, null)(TYPE(igraph_array3) *a);
IGRAPH_EXPORT BASE FUNCTION(igraph_array3, sum)(const TYPE(igraph_array3) *a);
IGRAPH_EXPORT void FUNCTION(igraph_array3, scale)(TYPE(igraph_array3) *a, BASE by);
IGRAPH_EXPORT void FUNCTION(igraph_array3, fill)(TYPE(igraph_array3) *a, BASE e);
IGRAPH_EXPORT int FUNCTION(igraph_array3, update)(TYPE(igraph_array3) *to,
                                                  const TYPE(igraph_array3) *from);
