/*
   IGraph library.
   Copyright (C) 2011-2012  Gabor Csardi <csardi@rmki.kfki.hu>
   334 Harvard street, Cambridge MA, 02139 USA

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

#include "igraph_vector.h"
#include "igraph_types_internal.h"

typedef struct {
    void *scanner;
    int eof;
    char errmsg[300];
    igraph_vector_t *vector;
    igraph_bool_t directed;
    int vcount, vcount2;
    int actfrom;
    int actto;
    int mode; /* 0: general, 1: vertex, 2: edge */
    igraph_trie_t *vertex_attribute_names;
    igraph_vector_ptr_t *vertex_attributes;
    igraph_trie_t *edge_attribute_names;
    igraph_vector_ptr_t *edge_attributes;
    int vertexid;
    int actvertex;
    int actedge;
} igraph_i_pajek_parsedata_t;
