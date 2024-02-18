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

#include "igraph_error.h"
#include "igraph_vector.h"
#include "igraph_vector_ptr.h"

#include "core/trie.h"

/* According to Pajek's author, limits of the Pajek program as of 2022-1-1 are:
 * "At the moment regular Pajek has limit one billion vertices,
 *  PajekXXL two billions, while Pajek 3XL ten billions."
 * Hard-coding the limit INT32_MAX is safe when compiling wiht 32-bit integers,
 * and likely sufficient for practical applications.
 */
#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
/* Limit maximum vertex count when using a fuzzer, to avoid out-of-memory failure. */
#define IGRAPH_PAJEK_MAX_VERTEX_COUNT (1 << 18)
#else
#define IGRAPH_PAJEK_MAX_VERTEX_COUNT INT32_MAX
#endif

#define CHECK_OOM_RP(p) IGRAPH_CHECK_OOM((p), "Not enough memory to read Pajek format.")
#define CHECK_OOM_WP(p) IGRAPH_CHECK_OOM((p), "Not enough memory to write Pajek format.")

typedef struct {
    void *scanner;
    igraph_bool_t eof;
    char errmsg[300];
    igraph_error_t igraph_errno;
    igraph_vector_int_t *vector;
    igraph_bool_t directed;
    igraph_integer_t vcount, vcount2;
    igraph_integer_t actfrom;
    igraph_integer_t actto;
    igraph_trie_t *vertex_attribute_names;
    igraph_vector_ptr_t *vertex_attributes;
    igraph_trie_t *edge_attribute_names;
    igraph_vector_ptr_t *edge_attributes;
    igraph_integer_t vertexid;
    igraph_integer_t actvertex;
    igraph_integer_t actedge;
} igraph_i_pajek_parsedata_t;
