/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_error.h"
#include "igraph_types.h"

#include "core/trie.h"

typedef enum { IGRAPH_DL_MATRIX,
               IGRAPH_DL_EDGELIST1, IGRAPH_DL_NODELIST1
             } igraph_i_dl_type_t;

typedef struct {
    void *scanner;
    int eof;
    int mode;
    long int n;
    long int from, to;
    igraph_vector_t edges;
    igraph_vector_t weights;
    igraph_strvector_t labels;
    igraph_trie_t trie;
    igraph_i_dl_type_t type;
    char errmsg[300];
} igraph_i_dl_parsedata_t;
