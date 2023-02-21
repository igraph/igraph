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

/* TODO: Find out maximum supported vertex count. */
#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
/* Limit maximum vertex count when using a fuzzer, to avoid out-of-memory failure. */
#define IGRAPH_DL_MAX_VERTEX_COUNT (1 << 20)
#else
#define IGRAPH_DL_MAX_VERTEX_COUNT INT32_MAX
#endif

typedef enum { IGRAPH_DL_MATRIX,
               IGRAPH_DL_EDGELIST1, IGRAPH_DL_NODELIST1
             } igraph_i_dl_type_t;

typedef struct {
    void *scanner;
    int eof;
    char errmsg[300];
    igraph_error_t igraph_errno;
    int mode;
    igraph_integer_t n;
    igraph_integer_t from, to;
    igraph_vector_int_t edges;
    igraph_vector_t weights;
    igraph_strvector_t labels;
    igraph_trie_t trie;
    igraph_i_dl_type_t type;
} igraph_i_dl_parsedata_t;
